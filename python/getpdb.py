import sys, requests, shutil
from modeller import *
from modeller.automodel import *

def pir(seq):
    s = "".join([l.strip() for l in seq.splitlines() if not l.startswith(">")])
    return f">P1;T\nsequence:T:::::::0.00:0.00\n{s}*"

def run(pdb):
    env = Environ()
    env.io.atom_files_directory = ["."]
    
    fasta = requests.get(f"https://www.rcsb.org/fasta/entry/{pdb}/download").text
    open("t.pir","w").write(pir(fasta))

    if not os.path.exists(pdb+".pdb"):
        open(pdb+".pdb","w").write(
            requests.get(f"https://files.rcsb.org/download/{pdb.upper()}.pdb").text
        )

    aln = Alignment(env)
    mdl = Model(env, file=pdb, model_segment=('FIRST:A','LAST:A'))
    aln.append_model(mdl, align_codes=pdb+"A")
    aln.append(file="t.pir", align_codes="T")

    aln.align2d()
    aln.write("a.ali")

    class M(AutoModel): pass
    m = M(env, alnfile="a.ali", knowns=pdb+"A", sequence="T")
    m.starting_model = m.ending_model = 1
    m.make()

    shutil.move("T.B99990001.pdb", "output.pdb")
    print("done")

if __name__ == "__main__":
    run(sys.argv[1])
