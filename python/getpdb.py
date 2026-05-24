import sys
import requests
from pymol import cmd
from modeller import *
from modeller.automodel import *
def build_missing_residues_model(pdbid: str, out_prefix: str):
    pdbid = pdbid.upper()
    print(f"Downloading PDB: {pdbid}")
    pdb_text = requests.get(
        f"https://files.rcsb.org/download/{pdbid}.pdb",
        timeout=60
    )
    pdb_text.raise_for_status()
    raw_pdb = f"{out_prefix}.pdb"
    with open(raw_pdb, "w") as f:
        f.write(pdb_text.text)
    print("Cleaning PDB with PyMOL...")
    cmd.load(raw_pdb, "model")
    cmd.remove("not polymer.protein")
    clean_pdb = f"{out_prefix}_clean.pdb"
    cmd.save(clean_pdb, "model")
    cmd.delete("all")
    print("Downloading FASTA...")
    fasta_text = requests.get(
        f"https://www.rcsb.org/fasta/entry/{pdbid}/download",
        timeout=60
    )
    fasta_text.raise_for_status()
    fasta_file = f"{out_prefix}.fasta"
    with open(fasta_file, "w") as f:
        f.write(fasta_text.text)
    print("Generating alignment...")
    log.verbose()
    env = Environ()
    env.io.atom_files_directory = ['.']
    aln = Alignment(env)
    mdl = Model(env, file=clean_pdb)
    aln.append_model(
        mdl,
        align_codes=pdbid,
        atom_files=clean_pdb
    )
    aln.append(
        file=fasta_file,
        align_codes="TARGET"
    )
    aln.align2d()
    aln_file = f"{out_prefix}.ali"
    aln.write(
        file=aln_file,
        alignment_format="PIR"
    )
    aln.write(
        file=f"{out_prefix}.pap",
        alignment_format="PAP"
    )
    print("Running MODELLER...")
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
    print("MODELLER finished")
def main():
    if len(sys.argv) != 3:
        print("Usage:")
        print("python script.py PDBID output_prefix")
        sys.exit(1)
    pdbid = sys.argv[1]
    out_prefix = sys.argv[2]
    build_missing_residues_model(
        pdbid,
        out_prefix
    )
if __name__ == "__main__":
    main()
