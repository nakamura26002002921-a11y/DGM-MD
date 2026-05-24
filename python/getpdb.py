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


def build_model(pdbid, out_prefix):

    pdbid = pdbid.upper()

    # -------------------------
    # 1. PDB download + clean
    # -------------------------
    pdb = requests.get(f"https://files.rcsb.org/download/{pdbid}.pdb")
    pdb.raise_for_status()

    raw_pdb = f"{out_prefix}.pdb"
    open(raw_pdb, "w").write(pdb.text)

    cmd.load(raw_pdb, "m")
    cmd.remove("not polymer.protein")
    clean_pdb = f"{out_prefix}_clean.pdb"
    cmd.save(clean_pdb, "m")
    cmd.delete("all")

    # -------------------------
    # 2. FASTA (target)
    # -------------------------
    fasta = requests.get(f"https://www.rcsb.org/fasta/entry/{pdbid}/download")
    fasta.raise_for_status()

    pir = fasta_to_pir(fasta.text, "TARGET")
    pir_file = f"{out_prefix}.pir"
    open(pir_file, "w").write(pir)

    # -------------------------
    # 3. MODELLER alignment
    # -------------------------
    env = Environ()
    env.io.atom_files_directory = ['.']
    env.io.hetatm = True   # iterative + ligand対応思想

    aln = Alignment(env)

    mdl = Model(env, file=clean_pdb, model_segment=('FIRST:A', 'LAST:A'))
    aln.append_model(mdl, align_codes=pdbid)

    aln.append(file=pir_file, align_codes="TARGET")

    # iterative tutorialの核心
    aln.align2d()

    aln_file = f"{out_prefix}.ali"
    aln.write(file=aln_file, alignment_format="PIR")

    # -------------------------
    # 4. Model building (iterative step 1)
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

    print("DONE: inspect DOPE and refine alignment if needed")


if __name__ == "__main__":
    build_model(sys.argv[1], sys.argv[2])
