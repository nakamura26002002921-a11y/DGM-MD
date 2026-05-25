import os,glob
from contextlib import contextmanager
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

@contextmanager
def _in_dir(p):
    o=os.getcwd();os.makedirs(p,exist_ok=True);os.chdir(p)
    try:yield p
    finally:os.chdir(o)

def fasta_to_pir(f,c):
    s=''.join(i.strip()for i in f.splitlines()if not i.startswith('>'))
    return f'>P1;{c}\nsequence:{c}:::::::0.00: 0.00\n{s}*\n'

def _best_dope_model(w,s):
    g=glob.glob(os.path.join(w,f'{s}.B9999*.pdb'))
    if not g:raise FileNotFoundError(f'モデル PDB が見つかりません: {os.path.join(w,f"{s}.B9999*.pdb")}')
    e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib');e.libs.parameters.read(file='$(LIB)/par.lib')
    b,k=None,float('inf')
    for p in g:
        d=Selection(complete_pdb(e,p)).assess_dope()
        if d<k:k,b=d,p
    return b

def build_profile(workdir,query_pir,pdb_db_pir,query_code='TARGET',matrix_offset=-450,gap_penalties_1d=(-500,-50),max_aln_evalue=.01,n_prof_iterations=1):
    with _in_dir(workdir):
        log.verbose();e=Environ();s=SequenceDB(e);b='pdb_95.bin'
        if not os.path.exists(b):
            s.read(seq_database_file=pdb_db_pir,seq_database_format='PIR',chains_list='ALL',minmax_db_seq_len=(30,4000),clean_sequences=True)
            s.write(seq_database_file=b,seq_database_format='BINARY',chains_list='ALL')
        s.read(seq_database_file=b,seq_database_format='BINARY',chains_list='ALL')
        a=Alignment(e);a.append(file=query_pir,alignment_format='PIR',align_codes='ALL')
        p=a.to_profile()
        p.build(s,matrix_offset=matrix_offset,rr_file='${LIB}/blosum62.sim.mat',gap_penalties_1d=gap_penalties_1d,n_prof_iterations=n_prof_iterations,check_profile=False,max_aln_evalue=max_aln_evalue)
        p.write(file='build_profile.prf',profile_format='TEXT')
        a=p.to_alignment();a.write(file='build_profile.ali',alignment_format='PIR')
    return os.path.join(workdir,'build_profile.prf')

def compare_templates(workdir,template_list,atom_files_dir='.'):
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[atom_files_dir,'.'];a=Alignment(e)
        for p,c in template_list:
            a.append_model(Model(e,file=p,model_segment=(f'FIRST:{c}',f'LAST:{c}')),atom_files=p,align_codes=p+c)
        a.malign();a.malign3d();a.compare_structures();a.id_table(matrix_file='family.mat');e.dendrogram(matrix_file='family.mat',cluster_cut=-1.)
    return os.path.join(workdir,'compare.log')

def align2d_single(workdir,template_pdb,template_chain,query_pir,query_code,template_code=None,max_gap_length=50):
    b=os.path.splitext(os.path.basename(template_pdb))[0]
    if template_code is None:template_code=b+template_chain
    o=f'{query_code}-{template_code}.ali'
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(template_pdb))or'.','.']
        a=Alignment(e)
        a.append_model(Model(e,file=template_pdb,model_segment=(f'FIRST:{template_chain}',f'LAST:{template_chain}')),align_codes=template_code,atom_files=template_pdb)
        a.append(file=query_pir,align_codes=query_code);a.align2d(max_gap_length=max_gap_length)
        a.write(file=o,alignment_format='PIR');a.write(file=o.replace('.ali','.pap'),alignment_format='PAP')
    return os.path.join(workdir,o)

def build_single_model(workdir,ali_file,template_code,sequence,n_models=5,assess_methods=None):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file))or'.','.']
        a=AutoModel(e,alnfile=os.path.abspath(ali_file),knowns=template_code,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    return sorted(glob.glob(os.path.join(workdir,f'{sequence}.B9999*.pdb')))

def evaluate_model(workdir,pdb_file,output_profile=None,smoothing_window=15):
    b=os.path.splitext(os.path.basename(pdb_file))[0]
    if output_profile is None:output_profile=f'{b}.profile'
    with _in_dir(workdir):
        log.verbose();e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib');e.libs.parameters.read(file='$(LIB)/par.lib')
        Selection(complete_pdb(e,os.path.abspath(pdb_file))).assess_dope(output='ENERGY_PROFILE NO_REPORT',file=output_profile,normalize_profile=True,smoothing_window=smoothing_window)
    return os.path.join(workdir,output_profile)

def salign_multiple_templates(workdir,template_list,atom_files_dir='.',output_ali='templates_mult.ali'):
    with _in_dir(workdir):
        log.verbose();e=Environ();e.io.atom_files_directory=[atom_files_dir,'.'];a=Alignment(e)
        for c,h in template_list:
            a.append_model(Model(e,file=c,model_segment=(f'FIRST:{h}',f'LAST:{h}')),atom_files=c,align_codes=c+h)
        for w,f,o in(((1.,0.,0.,0.,1.,0.),False,True),((1.,.5,1.,1.,1.,0.),False,True),((1.,1.,1.,1.,1.,0.),True,False)):
            a.salign(rms_cutoff=3.5,normalize_pp_scores=False,rr_file='$(LIB)/as1.sim.mat',overhang=30,gap_penalties_1d=(-450,-50),gap_penalties_3d=(0,3),gap_gap_score=0,gap_residue_score=0,dendrogram_file='templates.tree',alignment_type='tree',feature_weights=w,improve_alignment=True,fit=True,write_fit=f,write_whole_pdb=o,output='ALIGNMENT QUALITY')
        a.write(file=output_ali.replace('.ali','.pap'),alignment_format='PAP');a.write(file=output_ali,alignment_format='PIR')
        a.salign(rms_cutoff=1.,normalize_pp_scores=False,rr_file='$(LIB)/as1.sim.mat',overhang=30,gap_penalties_1d=(-450,-50),gap_penalties_3d=(0,3),gap_gap_score=0,gap_residue_score=0,dendrogram_file='templates_quality.tree',alignment_type='progressive',feature_weights=[0]*6,improve_alignment=False,fit=False,write_fit=True,write_whole_pdb=False,output='QUALITY')
    return os.path.join(workdir,output_ali)

def align2d_multiple(workdir,templates_ali,query_pir,query_code,output_ali='query_mult.ali',max_gap_length=20):
    with _in_dir(workdir):
        log.verbose();e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib')
        a=Alignment(e);a.append(file=os.path.abspath(templates_ali),align_codes='all');b=len(a)
        a.append(file=os.path.abspath(query_pir),align_codes=query_code)
        a.salign(output='',max_gap_length=max_gap_length,gap_function=True,alignment_type='PAIRWISE',align_block=b,feature_weights=(1.,0.,0.,0.,0.,0.),overhang=0,gap_penalties_1d=(-450,0),gap_penalties_2d=(.35,1.2,.9,1.2,.6,8.6,1.2,0.,0.),similarity_flag=True)
        a.write(file=output_ali,alignment_format='PIR');a.write(file=output_ali.replace('.ali','.pap'),alignment_format='PAP')
    return os.path.join(workdir,output_ali)

def build_multi_model(workdir,ali_file,template_codes,sequence,n_models=5,assess_methods=None,use_hetatm=False):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file))or'.','.'];e.io.hetatm=use_hetatm
        a=AutoModel(e,alnfile=os.path.abspath(ali_file),knowns=template_codes,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    return sorted(glob.glob(os.path.join(workdir,f'{sequence}.B9999*.pdb')))

def refine_loop(workdir,ini_model,sequence,loop_range,n_loop_models=10,md_level=None,atom_files_dir='.'):
    if md_level is None:md_level=refine.very_fast
    s,e=loop_range
    class _MyLoop(LoopModel):
        def select_loop_atoms(self):return Selection(self.residue_range(s,e))
    with _in_dir(workdir):
        log.verbose();v=Environ();v.io.atom_files_directory=[atom_files_dir,'.']
        m=_MyLoop(v,inimodel=os.path.abspath(ini_model),sequence=sequence)
        m.loop.starting_model=1;m.loop.ending_model=n_loop_models;m.loop.md_level=md_level;m.make()
    return sorted(glob.glob(os.path.join(workdir,f'{sequence}.BL*.pdb')))

def build_model_with_ligand(workdir,ali_file,template_codes,sequence,restraint_atom_pairs=None,restraint_mean=3.5,restraint_stdev=.1,n_models=5):
    if restraint_atom_pairs is None:restraint_atom_pairs=[]
    a,m,s=restraint_atom_pairs,restraint_mean,restraint_stdev
    class _MyModel(AutoModel):
        def special_restraints(self,aln):
            r=self.restraints
            for i,j in a:r.add(forms.UpperBound(group=physical.upper_distance,feature=features.Distance(self.atoms[i],self.atoms[j]),mean=m,stdev=s))
    with _in_dir(workdir):
        e=Environ();e.io.hetatm=True;e.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file))or'.','.']
        x=_MyModel(e,alnfile=os.path.abspath(ali_file),knowns=template_codes,sequence=sequence)
        x.starting_model=1;x.ending_model=n_models;x.make()
    return sorted(glob.glob(os.path.join(workdir,f'{sequence}.B9999*.pdb')))

def align2d_with_ss(workdir,template_pdb,template_chain,query_pir,query_code,output_ali=None,max_gap_length=50):
    b=os.path.splitext(os.path.basename(template_pdb))[0];t=b+template_chain
    if output_ali is None:output_ali=f'{query_code}-{t}.ali'
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(template_pdb))or'.','.']
        a=Alignment(e)
        a.append_model(Model(e,file=template_pdb,model_segment=(f'FIRST:{template_chain}',f'LAST:{template_chain}')),align_codes=t,atom_files=template_pdb)
        a.append(file=os.path.abspath(query_pir),align_codes=query_code);a.align2d(max_gap_length=max_gap_length)
        a.write(file=output_ali,alignment_format='PIR')
        a.write(file=output_ali.replace('.ali','.pap'),alignment_format='PAP',alignment_features='INDICES HELIX BETA')
    return os.path.join(workdir,output_ali)

def iterative_modeling(workdir,template_pdb,template_chain,query_pir,query_code,max_iterations=5,n_models=1):
    b=os.path.splitext(os.path.basename(template_pdb))[0];t=b+template_chain
    k,p=float('inf'),None
    for i in range(1,max_iterations+1):
        d=os.path.join(workdir,f'iter_{i:02d}');os.makedirs(d,exist_ok=True)
        print(f'\n=== 反復 {i}/{max_iterations} ===')
        a=align2d_with_ss(workdir=d,template_pdb=template_pdb,template_chain=template_chain,query_pir=query_pir,query_code=query_code)
        g=build_single_model(workdir=d,ali_file=a,template_code=t,sequence=query_code,n_models=n_models,assess_methods=(assess.DOPE,assess.GA341))
        if not g:
            print(f'  モデル生成に失敗しました (反復 {i})');break
        e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib');e.libs.parameters.read(file='$(LIB)/par.lib')
        q,r=float('inf'),None
        for x in g:
            y=Selection(complete_pdb(e,x)).assess_dope()
            if y<q:q,r=y,x
        print(f'  最良 DOPE スコア: {q:.2f}  ({r})')
        if q<k:k,p=q,r
        else:
            print('  スコアが改善しませんでした。反復を終了します。');break
    print(f'\n最終ベストモデル: {p} (DOPE={k:.2f})')
    return p

def build_from_threading_ali(workdir,ali_file,template_code,sequence,n_models=5,assess_methods=None):
    if assess_methods is None:assess_methods=(assess.DOPE,assess.GA341)
    with _in_dir(workdir):
        e=Environ();e.io.atom_files_directory=[os.path.dirname(os.path.abspath(ali_file))or'.','.']
        a=AutoModel(e,alnfile=os.path.abspath(ali_file),knowns=template_code,sequence=sequence,assess_methods=assess_methods)
        a.starting_model=1;a.ending_model=n_models;a.make()
    p=sorted(glob.glob(os.path.join(workdir,f'{sequence}.B9999*.pdb')))
    e=Environ();e.libs.topology.read(file='$(LIB)/top_heav.lib');e.libs.parameters.read(file='$(LIB)/par.lib')
    s=[]
    for i in p:s.append((Selection(complete_pdb(e,i)).assess_dope(),i))
    s.sort()
    return[i for _,i in s]

if __name__=="__main__":
    import sys,os,requests
    if len(sys.argv)<2:print(f"usage: python3 {sys.argv[0]} <PDB_ID>");sys.exit(1)
    t=sys.argv[1].upper();w=f"./modeller_work/{t}";os.makedirs(w,exist_ok=True)
    p=f"{w}/{t}.pdb";f=f"{w}/{t}.fasta";r=f"{w}/{t}.pir"
    if not os.path.exists(p):x=requests.get(f"https://files.rcsb.org/download/{t}.pdb",timeout=60);x.raise_for_status();open(p,"wb").write(x.content)
    if not os.path.exists(f):x=requests.get(f"https://www.rcsb.org/fasta/entry/{t}/download",timeout=60);x.raise_for_status();open(f,"w").write(x.text)
    if not os.path.exists(r):open(r,"w").write(fasta_to_pir(open(f).read(),t))
    a=align2d_single(workdir=w,template_pdb=p,template_chain='A',query_pir=r,query_code=t)
    build_single_model(workdir=w,ali_file=a,template_code=f"{t}A",sequence=t,n_models=5)
    print(_best_dope_model(w,t))
