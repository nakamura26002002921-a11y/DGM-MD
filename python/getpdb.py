import os, sys, requests
from modeller import *
from modeller.automodel import *

def fasta_to_pir(t, c):
    s = "".join([l.strip() for l in t.splitlines() if not l.startswith(">")])
    return f">P1;{c}\nsequence:{c}:::::::0.00:0.00\n{s}*"

def build_model(target_id, template_ids):
    env = Environ()
    env.io.atom_files_directory = ["."]
    
    # 1. ターゲット配列の取得
    fasta_url = f"https://www.rcsb.org/fasta/entry/{target_id}/download"
    fasta_text = requests.get(fasta_url).text
    with open("target.pir", "w") as f:
        f.write(fasta_to_pir(fasta_text, "TARGET"))

    # 2. テンプレートの比較と選定
    aln = Alignment(env)
    for pdb in template_ids:
        pdb = pdb.lower()
        # PDB取得（未ダウンロードの場合）
        if not os.path.exists(f"{pdb}.pdb"):
            r = requests.get(f"https://files.rcsb.org/download/{pdb.upper()}.pdb")
            with open(f"{pdb}.pdb", "w") as f: f.write(r.text)
        
        # モデルの読み込み (チェーンAを想定)
        m = Model(env, file=pdb, model_segment=('FIRST:A', 'LAST:A'))
        aln.append_model(m, atom_files=pdb, align_codes=pdb+"A")

    # 構造の重ね合わせと比較
    aln.malign()
    aln.malign3d()
    aln.compare_structures()
    aln.id_table(matrix_file='family.mat')
    env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)
    
    print("--- 構造比較が完了しました。family.mat とログを確認してください ---")
    
    # ここでは便宜上、最初のテンプレート（index 0）を選択してモデリングへ進みます
    best_template = template_ids[0].lower() + "A"
    
    # 3. 選定したテンプレートでのアライメントとモデリング
    aln = Alignment(env)
    mdl = Model(env, file=best_template[:-1], model_segment=('FIRST:A', 'LAST:A'))
    aln.append_model(mdl, atom_files=best_template[:-1], align_codes=best_template)
    aln.append(file="target.pir", align_codes="TARGET")
    
    aln.align2d(max_gap_length=50)
    aln.write(file="alignment.ali", alignment_format="PIR")
    
    class M(AutoModel): pass
    a = M(env, alnfile="alignment.ali", knowns=best_template, sequence="TARGET")
    a.starting_model = 1
    a.ending_model = 1
    a.make()
    
    print(f"DONE: Model built using {best_template}")

if __name__ == "__main__":
    # 例: python script.py 1ABC 1XYZ,2PDB,3GHI
    if len(sys.argv) != 3:
        print("Usage: python script.py <TARGET_ID> <TEMP1,TEMP2,...>")
        sys.exit(1)
    target = sys.argv[1]
    templates = sys.argv[2].split(",")
    build_model(target, templates)
