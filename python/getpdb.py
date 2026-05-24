import os,sys,requests
from pymol import cmd
from modeller import *
from modeller.automodel import *

def fasta_to_pir(t,c):
    s="".join([l.strip() for l in t.splitlines() if not l.startswith(">")])
    return f">P1;{c}\nsequence:{c}:::::::0.00:0.00\n{s}*"

def build_model(p):
    p=p.upper()
    u=f"https://files.rcsb.org/download/{p}.pdb"
    r=requests.get(u);r.raise_for_status()
    open("input.pdb","w").write(r.text)
    cmd.load("input.pdb","m")
    cmd.remove("not polymer.protein")
    cmd.save("template_clean.pdb","m")
    cmd.delete("all")
    f=f"https://www.rcsb.org/fasta/entry/{p}/download"
    r=requests.get(f);r.raise_for_status()
    open("target.pir","w").write(fasta_to_pir(r.text,"TARGET"))
    e=Environ();e.io.atom_files_directory=["."];e.io.hetatm=True
    a=Alignment(e)
    m=Model(e,file="template_clean.pdb",model_segment=('FIRST:A','LAST:A'))
    a.append_model(m,align_codes="TEMPLATE")
    a.append(file="target.pir",align_codes="TARGET")
    a.align2d()
    a.write(file="alignment.ali",alignment_format="PIR")
    class M(AutoModel): pass
    x=M(e,alnfile="alignment.ali",knowns="TEMPLATE",sequence="TARGET")
    x.starting_model=1;x.ending_model=1
    x.make()
    o="TARGET.B99990001.pdb"
    if os.path.exists(o):os.rename(o,"clean.pdb")
    print("DONE")

if __name__=="__main__":
    if len(sys.argv)!=2:sys.exit(1)
    build_model(sys.argv[1])
