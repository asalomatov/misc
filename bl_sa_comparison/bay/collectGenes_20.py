genes_20_blood = set()
genes_20_saliva = set()
import pandas as pd
import numpy as np
import os
from pyensembl import EnsemblRelease
data = EnsemblRelease(75)

df_descr = pd.read_csv("/mnt/xfs1/scratch/asalomatov/BloodVsSaliva/BloodVsSaliva_Descr.csv") 
df_descr.head()
for c in df_descr.columns:
    df_descr[c] = df_descr[c].map(str.strip)

df_descr['fam'] = df_descr['fam_id'].apply(lambda x: x.split("-")[0])
df_descr['protocol'] = df_descr['fam_id'].apply(lambda x: x.split("-")[1])

for smpl in df_descr['smpl_id']:
    print smpl
   # fam = df_descr['fam_id'][df_descr.smpl_id.isin([smpl])].iloc[0]
   # rel = df_descr['relationship'][df_descr.smpl_id.isin([smpl])].iloc[0]
    prot = df_descr['protocol'][df_descr.smpl_id.isin([smpl])].iloc[0]
    print prot
   # desc = '_'.join([fam, rel])
   # print desc
    fname='/mnt/xfs1/scratch/asalomatov/BloodVsSaliva/output/'+smpl+'-D1.final-exome-cvrg.txt'
    if not os.path.isfile(fname):
        print fname
        continue
    df_exome = pd.read_csv(fname, sep="\t", header=None) 
    df_exome = df_exome[df_exome[0] != 'all']
    df_exome.columns = ['chr', 'start', 'end', 'bin', 'length', 'ex_length', 'perc']
    df_exome['w'] = df_exome['perc']*df_exome['bin']
    df_exome['chr_start'] = df_exome['chr'].astype(str)+'_'+df_exome['start'].astype(str)
    df_ex_sum = df_exome.groupby(['chr_start']).w.sum()
    c1 = df_ex_sum < 20
    mygenes = set()
    for i, n in df_ex_sum[c1].iteritems():
        cont, sta = i.split("_")
        g_list = data.gene_names_at_locus(contig=cont, position=int(sta))
        for g in g_list:
            mygenes.add(g)
    if prot == 'blood':
        print 'adding to blood'
        genes_20_blood = genes_20_blood.union(mygenes)
    else:
        print 'adding to saliva'
        genes_20_saliva = genes_20_saliva.union(mygenes)
    print len(genes_20_blood), len(genes_20_saliva)


