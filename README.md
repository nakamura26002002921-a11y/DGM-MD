node.json をロード
edge.json をロード
while True
    現在状態を推定
        └ node の ls を順番に実行
    終了条件判定
    現在ノードから edge 選択
        └ weight 最小
    edge の command 実行
    出力ファイル確認
    if 成功
        retries reset
        次ループへ
    else
        retry count++
        if retry over
            abort
        MCTS repair
            sandbox 内で
                command 実行
                出力確認
            成功 command が見つかれば
                edge.json 更新
            else
                元 command 継続
