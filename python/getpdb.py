import os
import sys
import requests
from pymol import cmd
from modeller import *
from modeller.automodel import *


def fasta_to_pir(fasta_text, code):
    seq = "".join([l.strip() for l in fasta_text.splitlines() if not l.startswith(">")])
    return f""">P1;{code}
sequence:{code}:::::::0.00:0.00
{seq}*
"""


def build_model(pdbid):

    pdbid = pdbid.upper()

    # -------------------------
    # 1. PDB download + clean
    # -------------------------
    pdb_url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    pdb = requests.get(pdb_url)
    pdb.raise_for_status()

    raw_pdb = "input.pdb"
    with open(raw_pdb, "w") as f:
        f.write(pdb.text)

    cmd.load(raw_pdb, "m")
    cmd.remove("not polymer.protein")
    clean_pdb = "template_clean.pdb"
    cmd.save(clean_pdb, "m")
    cmd.delete("all")

    # -------------------------
    # 2. FASTA (target)
    # -------------------------
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdbid}/download"
    fasta = requests.get(fasta_url)
    fasta.raise_for_status()

    pir = fasta_to_pir(fasta.text, "TARGET")
    pir_file = "target.pir"
    with open(pir_file, "w") as f:
        f.write(pir)

    # -------------------------
    # 3. MODELLER setup
    # -------------------------
    env = Environ()
    env.io.atom_files_directory = ['.']
    env.io.hetatm = True

    aln = Alignment(env)

    mdl = Model(env, file=clean_pdb, model_segment=('FIRST:A', 'LAST:A'))
    aln.append_model(mdl, align_codes=pdbid)

    aln.append(file=pir_file, align_codes="TARGET")

    aln.align2d()

    aln_file = "alignment.ali"
    aln.write(file=aln_file, alignment_format="PIR")

    # -------------------------
    # 4. Build model
    # -------------------------
    class MyModel(AutoModel):
        pass

    a = MyModel(
        env,
        alnfile=aln_file,
        knowns=pdbid,
        sequence="TARGET"
    )

    a.starting_model = 1
    a.ending_model = 1
    a.make()

    # -------------------------
    # 5. FINAL OUTPUT -> clean.pdb
    # -------------------------
    out_pdb = f"{pdbid}.B99990001.pdb"

    if os.path.exists(out_pdb):
        os.rename(out_pdb, "clean.pdb")

    print("DONE -> clean.pdb generated")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 main.py PDBID")
        sys.exit(1)

    build_model(sys.argv[1])
