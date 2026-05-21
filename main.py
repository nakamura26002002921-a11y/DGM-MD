#!/usr/bin/env python3
"""
MD Simulation Workflow Automation
 - 有向グラフに沿ってコマンドを逐次実行
 - 出力ファイルが生成されなかった場合、LLM支援のMCTSでコマンドを修復
"""

import json
import os
import re
import glob
import math
import random
import shutil
import logging
import tempfile
import subprocess
from copy import deepcopy
from datetime import datetime

# ── Anthropic (オプション) ──────────────────────────────────────────────────
try:
    import anthropic
    _client = anthropic.Anthropic()
    HAS_LLM = True
except Exception:
    _client = None
    HAS_LLM = False

# ── ディレクトリ準備 ─────────────────────────────────────────────────────────
LOG_DIR    = "logs"
SANDBOX_DIR = "sandbox"
os.makedirs(LOG_DIR,     exist_ok=True)
os.makedirs(SANDBOX_DIR, exist_ok=True)
os.makedirs("sys",       exist_ok=True)

_ts = datetime.now().strftime("%Y%m%d_%H%M%S")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"{LOG_DIR}/main_{_ts}.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

# ── ユーザー設定 ─────────────────────────────────────────────────────────────
CONFIG: dict[str, str] = {
    "PDBID":        os.environ.get("PDBID",  "1AKI"),
    "GMX":          os.environ.get("GMX",    "gmx"),
    "PYMOL":        os.environ.get("PYMOL",  "pymol"),
    "cle_pdb":      "sys/clean.pdb",
    "pro_gro":      "sys/protein.gro",
    "pro_top":      "sys/topol.top",
    "WATER_MODEL":  "tip3p",
    "FF":           "amber03",
    "DISTANCE":     "1.0",
    "box_gro":      "sys/newbox.gro",
    "WATERBOXFILE": "spc216.gro",
    "sol_gro":      "sys/solv.gro",
    "sol_top":      "sys/topol.top",
    "ions_mdp":     "ions.mdp",
    "ions_tpr":     "sys/ions.tpr",
    "ions_gro":     "sys/ions.gro",
    "ions_top":     "sys/topol.top",
    "tmp_pdb":      "sys/tmp.pdb",
}

TERMINAL_NODE   = 6   # sys/MD.gro が存在すれば完了
MCTS_ITERATIONS = 50  # MCTSの反復回数
MAX_RETRIES     = 3   # エッジごとの最大リトライ数


# ═══════════════════════════════════════════════════════════════
# グラフ I/O
# ═══════════════════════════════════════════════════════════════

def load_graph() -> tuple[dict, dict]:
    with open("node.json") as f:
        raw_nodes = json.load(f)
    with open("edge.json") as f:
        raw_edges = json.load(f)
    nodes = {int(k): v for k, v in raw_nodes.items()}
    edges = {int(k): v for k, v in raw_edges.items()}
    return nodes, edges


def save_edges(edges: dict) -> None:
    serialisable = {str(k): v for k, v in edges.items()}
    with open("edge.json", "w") as f:
        json.dump(serialisable, f, indent=4, ensure_ascii=False)
    log.info("edge.json を更新しました。")


# ═══════════════════════════════════════════════════════════════
# テンプレート展開・コマンド実行
# ═══════════════════════════════════════════════════════════════

def render_cmd(template: str, cfg: dict = CONFIG) -> str:
    """CONFIG の値で {KEY} プレースホルダを置換する。"""
    cmd = template
    for k, v in cfg.items():
        cmd = cmd.replace("{" + k + "}", str(v))
    return cmd


def run_cmd(
    cmd: str,
    label: str = "",
    cwd: str = ".",
) -> tuple[int, str, str]:
    """
    シェルコマンドを実行し (returncode, stdout, stderr) を返す。
    label が指定された場合はログファイルも保存する。
    """
    log.info(f"[RUN] {cmd}")
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, cwd=cwd
    )
    combined = (
        f"CMD: {cmd}\n"
        f"=== STDOUT ===\n{result.stdout}\n"
        f"=== STDERR ===\n{result.stderr}"
    )
    if label:
        log_path = os.path.join(LOG_DIR, f"{label}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        with open(log_path, "w") as lf:
            lf.write(combined)

    if result.returncode != 0:
        log.warning(f"[FAIL rc={result.returncode}] {cmd}")
    return result.returncode, result.stdout, result.stderr


# ═══════════════════════════════════════════════════════════════
# 状態推定
# ═══════════════════════════════════════════════════════════════

def estimate_state(nodes: dict) -> int:
    """
    各ノードの ls チェックを実行し、
    パスした最大のノードインデックスを現在の状態とする。
    """
    current = -1
    for idx in sorted(nodes.keys()):
        rc, _, _ = run_cmd(nodes[idx]["cmd"], label=f"state_node{idx}")
        if rc == 0:
            current = idx
    log.info(f"現在の状態: node {current}")
    return current


# ═══════════════════════════════════════════════════════════════
# エッジ選択・出力確認
# ═══════════════════════════════════════════════════════════════

def select_edge(current_node: int, edges: dict) -> int | None:
    """
    current_node を起点とするエッジの中から
    weight が最小のものを選択して返す。
    """
    candidates = [
        eid for eid, e in edges.items() if e["in"] == current_node
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda eid: edges[eid].get("weight", 1))


def output_exists(out_node: int, nodes: dict) -> bool:
    """
    出力ノードの ls コマンドが成功すれば True。
    ノードが未定義の場合は True（スキップ）とする。
    """
    if out_node not in nodes:
        log.debug(f"node {out_node} は node.json に未定義 → 出力チェックをスキップ")
        return True
    rc, _, _ = run_cmd(nodes[out_node]["cmd"], label=f"out_check_node{out_node}")
    return rc == 0


# ═══════════════════════════════════════════════════════════════
# LLM ヘルパー
# ═══════════════════════════════════════════════════════════════

def _collect_recent_logs(n: int = 2) -> str:
    """直近 n 個のログファイルの内容を返す。"""
    logs = sorted(
        glob.glob(f"{LOG_DIR}/*.log"),
        key=os.path.getmtime,
        reverse=True,
    )
    parts = []
    for lp in logs[:n]:
        try:
            with open(lp) as f:
                parts.append(f"=== {os.path.basename(lp)} ===\n{f.read()}")
        except Exception:
            pass
    return "\n\n".join(parts)


def llm_suggest_variants(failed_cmd: str, n: int = 5) -> list[str]:
    """
    失敗したコマンドとログを Claude に渡し、
    修正候補のコマンドリストを取得する。
    """
    if not HAS_LLM:
        return []
    recent_logs = _collect_recent_logs(2)
    prompt = (
        f"以下のシェルコマンドが GROMACS MD セットアップパイプラインで失敗しました:\n\n"
        f"```\n{failed_cmd}\n```\n\n"
        f"直近のログ:\n```\n{recent_logs[-3000:]}\n```\n\n"
        f"エラーを修正できる可能性のあるコマンドを {n} 個提案してください。\n"
        f"JSON 配列のみ返してください（前置き・コードフェンス不要）。\n"
        f"例: [\"cmd1\", \"cmd2\"]"
    )
    try:
        resp = _client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=1000,
            messages=[{"role": "user", "content": prompt}],
        )
        text = resp.content[0].text.strip()
        text = re.sub(r"^```[^\n]*\n?", "", text)
        text = re.sub(r"\n?```$", "", text.strip())
        variants = json.loads(text)
        if isinstance(variants, list):
            return [str(v) for v in variants[:n]]
    except Exception as exc:
        log.warning(f"LLM 提案失敗: {exc}")
    return []


# ═══════════════════════════════════════════════════════════════
# MCTS
# ═══════════════════════════════════════════════════════════════

class _MCTSNode:
    """MCTS の1ノード = 1つのコマンド候補。"""

    def __init__(self, cmd: str, parent: "_MCTSNode | None" = None):
        self.cmd      = cmd
        self.parent   = parent
        self.children: list["_MCTSNode"] = []
        self.visits   = 0
        self.wins     = 0

    def uct(self, c: float = 1.41) -> float:
        if self.visits == 0:
            return float("inf")
        exploit = self.wins / self.visits
        explore = c * math.sqrt(math.log(self.parent.visits) / self.visits)
        return exploit + explore

    def best_child(self) -> "_MCTSNode":
        return max(self.children, key=lambda n: n.uct())

    def is_leaf(self) -> bool:
        return len(self.children) == 0


def _sandbox_run(cmd: str, out_node: int, nodes: dict) -> bool:
    """
    一時ディレクトリ内でコマンドを試し、
    出力ノードの ls チェックが通れば True を返す。
    """
    sbox = tempfile.mkdtemp(dir=SANDBOX_DIR)
    try:
        # sys/ ディレクトリを丸ごとコピー
        sys_dst = os.path.join(sbox, "sys")
        if os.path.isdir("sys"):
            shutil.copytree("sys", sys_dst)
        else:
            os.makedirs(sys_dst)

        # mdp / itp ファイルもコピー
        for pattern in ("*.mdp", "*.itp", "*.top"):
            for fp in glob.glob(pattern):
                shutil.copy2(fp, sbox)

        rc, _, _ = run_cmd(cmd, cwd=sbox)
        if rc != 0:
            return False

        # 出力ノードのチェック（サンドボックス内で評価）
        if out_node in nodes:
            ls_cmd = nodes[out_node]["cmd"]
            rc2, _, _ = run_cmd(ls_cmd, cwd=sbox)
            return rc2 == 0
        return True

    finally:
        shutil.rmtree(sbox, ignore_errors=True)


def mcts_fix(
    failed_cmd: str,
    edge_idx: int,
    edges: dict,
    nodes: dict,
    out_node: int,
    n_iter: int = MCTS_ITERATIONS,
) -> str | None:
    """
    サンドボックス内でモンテカルロ木探索を実行し、
    成功したコマンド文字列を返す（見つからなければ None）。

    フロー:
      選択 → 展開 → シミュレーション → バックプロパゲーション
    """
    log.info(f"[MCTS] edge {edge_idx} の修復を開始 (iterations={n_iter}) …")

    # 初期候補: 失敗したコマンド + LLM 提案
    initial_variants = [failed_cmd] + llm_suggest_variants(failed_cmd, n=5)

    root = _MCTSNode(cmd=failed_cmd)
    root.visits   = 1
    root.children = [_MCTSNode(cmd=v, parent=root) for v in initial_variants]

    best_cmd: str | None = None

    for iteration in range(n_iter):
        # ── 選択 ─────────────────────────────────────────────────────────
        node = root
        while not node.is_leaf():
            node = node.best_child()

        # ── 展開 ─────────────────────────────────────────────────────────
        if node.visits > 2 and not node.children:
            new_variants = llm_suggest_variants(node.cmd, n=3)
            node.children = [_MCTSNode(cmd=v, parent=node) for v in new_variants]
            if node.children:
                node = random.choice(node.children)

        # ── シミュレーション ───────────────────────────────────────────
        success = _sandbox_run(node.cmd, out_node, nodes)

        # ── バックプロパゲーション ─────────────────────────────────────
        cur = node
        while cur is not None:
            cur.visits += 1
            if success:
                cur.wins += 1
            cur = cur.parent

        if success:
            log.info(f"[MCTS] iter={iteration} SUCCESS → {node.cmd}")
            best_cmd = node.cmd
            break
        else:
            log.debug(f"[MCTS] iter={iteration} fail")

    if best_cmd is None:
        log.warning("[MCTS] 有効なコマンドが見つかりませんでした。")

    return best_cmd


# ═══════════════════════════════════════════════════════════════
# メインループ
# ═══════════════════════════════════════════════════════════════

def main() -> None:
    log.info("=" * 60)
    log.info("MD ワークフロー 開始")
    log.info("=" * 60)

    nodes, edges = load_graph()
    retries: dict[int, int] = {}   # edge_id → 連続失敗回数

    while True:

        # ── 1. 現在の状態を推定 ───────────────────────────────────────
        current = estimate_state(nodes)

        # ── 2. 終了条件 ───────────────────────────────────────────────
        if current >= TERMINAL_NODE:
            log.info(f"終了ノード {TERMINAL_NODE} に到達。ワークフロー完了！")
            break

        # ── 3. エッジ選択 ─────────────────────────────────────────────
        edge_id = select_edge(current, edges)
        if edge_id is None:
            log.error(f"node {current} からの出力エッジがありません。停止します。")
            break

        edge     = edges[edge_id]
        out_node = edge["out"]
        log.info(f"node {current} → edge {edge_id} → node {out_node}")

        # ── 4. コマンド実行 ───────────────────────────────────────────
        cmd = render_cmd(edge["cmd"])
        run_cmd(cmd, label=f"edge{edge_id}")

        # ── 5. 出力ファイルの確認 ─────────────────────────────────────
        if output_exists(out_node, nodes):
            log.info(f"node {out_node} の出力を確認。状態を更新します。")
            retries[edge_id] = 0
            continue   # ループ先頭で状態を再推定

        # ── 6. 失敗処理 ───────────────────────────────────────────────
        retries[edge_id] = retries.get(edge_id, 0) + 1
        log.warning(
            f"node {out_node} の出力が確認できません。"
            f" (リトライ {retries[edge_id]}/{MAX_RETRIES})"
        )

        if retries[edge_id] > MAX_RETRIES:
            log.error(f"edge {edge_id} が最大リトライ数を超えました。中断します。")
            break

        # ── 7. MCTS による修復 ────────────────────────────────────────
        fixed_cmd = mcts_fix(cmd, edge_id, edges, nodes, out_node)

        if fixed_cmd:
            log.info(f"[MCTS] edge {edge_id} のコマンドを更新: {fixed_cmd}")
            edges[edge_id] = deepcopy(edge)
            edges[edge_id]["cmd"] = fixed_cmd
            save_edges(edges)
            retries[edge_id] = 0
        else:
            log.warning("[MCTS] 修復コマンドが見つかりませんでした。元のコマンドで再試行します。")

    log.info("=" * 60)
    log.info("MD ワークフロー 終了")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
