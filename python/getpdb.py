import os,glob
from contextlib import contextmanager
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

@contextmanager
def _in_dir(path):
    o=os.getcwd();os.makedirs(path,exist_ok=True);os.chdir(path)
    try:yield path
    finally:os.chdir(o)

def fasta_to_pir(fasta_text,code):
    seq="".join(l.strip() for l in fasta_text.splitlines() if not l.startswith(">"))
    return f">P1;{code}\nsequence:{code}:::::::0.00: 0.00\n{seq}*\n"

def _best_dope_model(workdir,sequence):
    pattern=os.path.join(workdir,f"{sequence}.B9999*.pdb")
    cs=glob.glob(pattern)
    if not cs:raise FileNotFoundError(pattern)
    env=Environ()
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    best,b=float("inf"),None
    for p in cs:
        mdl=complete_pdb(env,p);s=Selection(mdl);sc=s.assess_dope()
        if sc<best:best,bb=sc,p
    return bb

def build_profile(workdir,query_pir,pdb_db_pir,query_code="TARGET",matrix_offset=-450,gap_penalties_1d=(-500,-50),max_aln_evalue=0.01,n_prof_iterations=1):
    with _in_dir(workdir):
        log.verbose();env=Environ()
        sdb=SequenceDB(env);bin_db="pdb_95.bin"
        if not os.path.exists(bin_db):
            sdb.read(seq_database_file=pdb_db_pir,seq_database_format="PIR",chains_list="ALL",minmax_db_seq_len=(30,4000),clean_sequences=True)
            sdb.write(seq_database_file=bin_db,seq_database_format="BINARY",chains_list="ALL")
        sdb.read(seq_database_file=bin_db,seq_database_format="BINARY",chains_list="ALL")
        aln=Alignment(env);aln.append(file=query_pir,alignment_format="PIR",align_codes="ALL")
        prf=aln.to_profile()
        prf.build(sdb,matrix_offset=matrix_offset,rr_file="${LIB}/blosum62.sim.mat",gap_penalties_1d=gap_penalties_1d,n_prof_iterations=n_prof_iterations,check_profile=False,max_aln_evalue=max_aln_evalue)
        prf.write(file="build_profile.prf",profile_format="TEXT")
    return os.path.join(workdir,"build_profile.prf")

def compare_templates(workdir,template_list,atom_files_dir="."):
    with _in_dir(workdir):
        env=Environ();env.io.atom_files_directory=[atom_files_dir,"."]
        aln=Alignment(env)
        for p,c in template_list:
            m=Model(env,file=p,model_segment=(f"FIRST:{c}",f"LAST:{c}"))
            aln.append_model(m,atom_files=p,align_codes=p+c)
        aln.malign();aln.malign3d();aln.compare_structures()
        aln.id_table(matrix_file="family.mat");env.dendrogram(matrix_file="family.mat",cluster_cut=-1.0)
    return os.path.join(workdir,"compare.log")

def align2d_single(workdir,template_pdb,template_chain,query_pir,query_code,template_code=None,max_gap_length=50):
    pdb_base=os.path.splitext(os.path.basename(template_pdb))[0]
    template_code=template_code or pdb_base+template_chain
    ali_out=f"{query_code}-{template_code}.ali"
    with _in_dir(workdir):
        env=Environ()
        env.io.atom_files_directory=[os.path.dirname(os.path.abspath(template_pdb)) or ".", "."]
        aln=Alignment(env)
        mdl=Model(env,file=template_pdb,model_segment=(f"FIRST:{template_chain}",f"LAST:{template_chain}"))
        aln.append_model(mdl,align_codes=template_code,atom_files=template_pdb)
        aln.append(file=query_pir,align_codes=query_code)
        aln.align2d(max_gap_length=max_gap_length)
        aln.write(file=ali_out,alignment_format="PIR")
    return os.path.join(workdir,ali_out)

def build_single_model(workdir,ali_file,template_code,sequence,n_models=5,assess_methods=None):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        env=Environ()
        env.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file)) or ".", "."]
        a=AutoModel(env,alnfile=os.path.abspath(ali_file),knowns=template_code,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    return sorted(glob.glob(os.path.join(workdir,f"{sequence}.B9999*.pdb")))

def evaluate_model(workdir,pdb_file,output_profile=None,smoothing_window=15):
    pdb_base=os.path.splitext(os.path.basename(pdb_file))[0]
    output_profile=output_profile or f"{pdb_base}.profile"
    with _in_dir(workdir):
        log.verbose();env=Environ()
        env.libs.topology.read(file="$(LIB)/top_heav.lib")
        env.libs.parameters.read(file="$(LIB)/par.lib")
        mdl=complete_pdb(env,os.path.abspath(pdb_file))
        s=Selection(mdl)
        s.assess_dope(output="ENERGY_PROFILE NO_REPORT",file=output_profile,normalize_profile=True,smoothing_window=smoothing_window)
    return os.path.join(workdir,output_profile)

def salign_multiple_templates(workdir,template_list,atom_files_dir=".",output_ali="templates_mult.ali"):
    with _in_dir(workdir):
        env=Environ();env.io.atom_files_directory=[atom_files_dir,"."]
        aln=Alignment(env)
        for c,ch in template_list:
            mdl=Model(env,file=c,model_segment=(f"FIRST:{ch}",f"LAST:{ch}"))
            aln.append_model(mdl,atom_files=c,align_codes=c+ch)
        for w,rf,wh in(((1.,0.,0.,0.,1.,0.),False,True),((1.,0.5,1.,1.,1.,0.),False,True),((1.,1.,1.,1.,1.,0.),True,False)):
            aln.salign(rms_cutoff=3.5,rr_file="$(LIB)/as1.sim.mat",overhang=30,gap_penalties_1d=(-450,-50),gap_penalties_3d=(0,3),feature_weights=w,improve_alignment=True,fit=True,write_fit=rf,write_whole_pdb=wh,output="ALIGNMENT QUALITY")
        aln.write(file=output_ali,alignment_format="PIR")
    return os.path.join(workdir,output_ali)

def align2d_multiple(workdir,templates_ali,query_pir,query_code,output_ali="query_mult.ali",max_gap_length=20):
    with _in_dir(workdir):
        env=Environ()
        aln=Alignment(env)
        aln.append(file=os.path.abspath(templates_ali),align_codes="all")
        n=len(aln)
        aln.append(file=os.path.abspath(query_pir),align_codes=query_code)
        aln.salign(output="",max_gap_length=max_gap_length,gap_function=True,alignment_type="PAIRWISE",align_block=n)
        aln.write(file=output_ali,alignment_format="PIR")
    return os.path.join(workdir,output_ali)

def build_multi_model(workdir,ali_file,template_codes,sequence,n_models=5,assess_methods=None,use_hetatm=False):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        env=Environ()
        env.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file)) or ".", "."]
        env.io.hetatm=use_hetatm
        a=AutoModel(env,alnfile=os.path.abspath(ali_file),knowns=template_codes,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    return sorted(glob.glob(os.path.join(workdir,f"{sequence}.B9999*.pdb")))

def refine_loop(workdir,ini_model,sequence,loop_range,n_loop_models=10,md_level=None,atom_files_dir="."):
    md_level=md_level or refine.very_fast
    s,e=loop_range
    class L(LoopModel):
        def select_loop_atoms(self):return Selection(self.residue_range(s,e))
    with _in_dir(workdir):
        env=Environ();env.io.atom_files_directory=[atom_files_dir,"."]
        m=L(env,inimodel=os.path.abspath(ini_model),sequence=sequence)
        m.loop.starting_model=1;m.loop.ending_model=n_loop_models;m.loop.md_level=md_level;m.make()
    return sorted(glob.glob(os.path.join(workdir,f"{sequence}.BL*.pdb")))

def build_model_with_ligand(workdir,ali_file,template_codes,sequence,restraint_atom_pairs=None,restraint_mean=3.5,restraint_stdev=0.1,n_models=5):
    restraint_atom_pairs=restraint_atom_pairs or []
    class M(AutoModel):
        def special_restraints(self,aln):
            rsr=self.restraints
            for i1,i2 in restraint_atom_pairs:
                a=[self.atoms[i1],self.atoms[i2]]
                rsr.add(forms.UpperBound(group=physical.upper_distance,feature=features.Distance(*a),mean=restraint_mean,stdev=restraint_stdev))
    with _in_dir(workdir):
        env=Environ();env.io.hetatm=True
        env.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file)) or ".", "."]
        a=M(env,alnfile=os.path.abspath(ali_file),knowns=template_codes,sequence=sequence)
        a.starting_model=1;a.ending_model=n_models;a.make()
    return sorted(glob.glob(os.path.join(workdir,f"{sequence}.B9999*.pdb")))

def align2d_with_ss(workdir,template_pdb,template_chain,query_pir,query_code,output_ali=None,max_gap_length=50):
    pdb_base=os.path.splitext(os.path.basename(template_pdb))[0]
    tc=pdb_base+template_chain
    output_ali=output_ali or f"{query_code}-{tc}.ali"
    with _in_dir(workdir):
        env=Environ()
        env.io.atom_files_directory=[os.path.dirname(os.path.abspath(template_pdb)) or ".", "."]
        mdl=Model(env,file=template_pdb,model_segment=(f"FIRST:{template_chain}",f"LAST:{template_chain}"))
        aln=Alignment(env)
        aln.append_model(mdl,align_codes=tc,atom_files=template_pdb)
        aln.append(file=os.path.abspath(query_pir),align_codes=query_code)
        aln.align2d(max_gap_length=max_gap_length)
        aln.write(file=output_ali,alignment_format="PIR")
    return os.path.join(workdir,output_ali)

def iterative_modeling(workdir,template_pdb,template_chain,query_pir,query_code,max_iterations=5,n_models=1):
    tc=os.path.splitext(os.path.basename(template_pdb))[0]+template_chain
    best=float("inf");bp=None
    for i in range(1,max_iterations+1):
        d=os.path.join(workdir,f"iter_{i:02d}");os.makedirs(d,exist_ok=True)
        ali=align2d_with_ss(d,template_pdb,template_chain,query_pir,query_code)
        pdbs=build_single_model(d,ali,tc,query_code,n_models)
        env=Environ();env.libs.topology.read(file="$(LIB)/top_heav.lib");env.libs.parameters.read(file="$(LIB)/par.lib")
        ib=float("inf");bp2=None
        for p in pdbs:
            m=complete_pdb(env,p);s=Selection(m);sc=s.assess_dope()
            if sc<ib:ib,bp2=sc,p
        if ib<best:best,bp=ib,bp2
        else:break
    return bp

def build_from_threading_ali(workdir,ali_file,template_code,sequence,n_models=5,assess_methods=None):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        env=Environ()
        env.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file)) or ".", "."]
        a=AutoModel(env,alnfile=os.path.abspath(ali_file),knowns=template_code,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    pdbs=sorted(glob.glob(os.path.join(workdir,f"{sequence}.B9999*.pdb")))
    env=Environ();env.libs.topology.read(file="$(LIB)/top_heav.lib");env.libs.parameters.read(file="$(LIB)/par.lib")
    sc=[]
    for p in pdbs:
        m=complete_pdb(env,p);s=Selection(m);sc.append((s.assess_dope(),p))
    sc.sort()
    return [p for _,p in sc]

if __name__=="__main__":
    import sys,os,requests
    t=sys.argv[1]
    w=f"./modeller_work/{t}"
    os.makedirs(w,exist_ok=True)
    p=f"{w}/{t}.pdb"
    r=f"{w}/{t}.pir"
    if not os.path.exists(p):open(p,"wb").write(requests.get(f"https://files.rcsb.org/download/{t}.pdb",timeout=60).content)
    if not os.path.exists(r):
        f=requests.get(f"https://www.rcsb.org/fasta/entry/{t}/download",timeout=60);f.raise_for_status()
        s="".join(i.strip() for i in f.text.splitlines() if not i.startswith(">"))
        open(r,"w").write(f">P1;{t}\nsequence:{t}:::::::0.00:0.00\n{s}*\n")
    ali=align2d_single(w,f"{t}.pdb","A",f"{t}.pir",t)
    pdbs=build_single_model(w,os.path.basename(ali),f"{t}A",t,5)
    if pdbs:open(f"{t}.pdb","wb").write(open(pdbs[0],"rb").read())
