import os,sys,requests
from pymol import cmd
from modeller import *
from modeller.automodel import *

def fasta_to_pir(t,c):
    s="".join([l.strip() for l in t.splitlines() if not l.startswith(">")])
    return f">P1;{c}\nsequence:{c}:::::::0.00:0.00\n{s}*"

import os
from modeller import *
from modeller.automodel import *

def run_modeller_pipeline(base_path, pdbid, target_code, sequence):
    pdbid = pdbid.upper()
    sys_dir = os.path.join(base_path, "sys")
    os.makedirs(sys_dir, exist_ok=True)
    
    # ファイルパス設定
    pir_file = os.path.join(sys_dir, "target.ali")
    ali_file = os.path.join(sys_dir, "alignment.ali")
    
    # 1. PIR形式の作成 (チュートリアル同様のヘッダーフォーマット)
    with open(pir_file, "w") as f:
        f.write(f">P1;{target_code}\nsequence:{target_code}:::::::0.00:0.00\n{sequence}*")

    env = Environ()
    env.io.atom_files_directory = [sys_dir]
    aln = Alignment(env)
    mdl = Model(env, file=pdbid)
    aln.append_model(mdl, align_codes=pdbid, atom_files=f"{pdbid}.pdb")
    aln.append(file=pir_file, align_codes=target_code)
    
    aln.align2d()
    aln.write(file=ali_file, alignment_format="PIR")
    
    # 4. モデル構築
    class M(AutoModel):
        def select_atoms(self):
            return Selection(self.chains)

    a = M(env, 
          alnfile=ali_file, 
          knowns=pdbid, 
          sequence=target_code)
    
    a.starting_model = 1
    a.ending_model = 1
    a.make()
    
    return f"{target_code}.B99990001.pdb"

if __name__=="__main__":
    if len(sys.argv)!=2:sys.exit(1)
    build_model(sys.argv[1])
