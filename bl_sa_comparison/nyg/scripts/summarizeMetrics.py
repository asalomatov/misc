j
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pylab
#pylab.show()
pd.__version__

### dict of sample - family
df_descr = pd.read_csv("/mnt/xfs1/scratch/asalomatov/BloodVsSaliva/BloodVsSaliva_Descr.csv") 
df_descr.head()
for c in df_descr.columns:
    df_descr[c] = df_descr[c].map(str.strip)

df_descr['fam'] = df_descr['fam_id'].apply(lambda x: x.split("-")[0])
df_descr['protocol'] = df_descr['fam_id'].apply(lambda x: x.split("-")[1])
df_descr.columns
df_descr['gender'][0]
df_descr.gender.value_counts()
df_descr.protocol.value_counts()

### read read counts
df_FlSt = pd.read_csv("/mnt/xfs1/scratch/asalomatov/BloodVsSaliva/output/metrics/BlVsSlv_Fl_total_reads.txt", sep=" ",
        header=None) 
df_FlSt.columns = ['f', 'subgroup', 'num_alignments', 'x1','x2']
df_FlSt.head()
df_FlSt.columns
df_FlSt['smpl_id'] = df_FlSt['f'].apply(lambda x: "-".join(x.split("-")[:-1]))
df_FlSt.head()
df_FlSt.tail()
df_FlSt.drop(['f', 'x1', 'x2'], axis=1, inplace=True)
df_FlSt = df_FlSt.merge(df_descr, 'left', on='smpl_id')
df_FlSt = df_FlSt.dropna()
df_FlSt['subgroup'][df_FlSt.subgroup.isin(['final'])] = 'total'
df_FlSt.protocol.value_counts()

### read mapped
df_mapped = pd.read_csv("/mnt/xfs1/scratch/asalomatov/BloodVsSaliva/output/metrics/BlVsSlv_Fl_mapped_reads.txt", sep=" ",
        header=None) 
df_mapped.columns = ['f', 'subgroup', 'num_mapped', 'x1','x2']
df_mapped.head()

df_mapped['smpl_id'] = df_mapped['f'].apply(lambda x: "-".join(x.split("-")[:-1]))
df_mapped.head()
df_mapped.tail()
df_mapped.drop(['f', 'x1', 'x2'], axis=1, inplace=True)
df_mapped = df_mapped.merge(df_descr, 'left', on='smpl_id')
df_mapped = df_mapped.dropna()
df_mapped['subgroup'][df_mapped.subgroup.isin(['final'])] = 'total'
df_mapped.protocol.value_counts()

### read duplicate
df_duplicate = pd.read_csv("/mnt/xfs1/scratch/asalomatov/BloodVsSaliva/output/metrics/BlVsSlv_Fl_duplicate_reads.txt", sep=" ",
        header=None) 
df_duplicate.columns = ['f', 'subgroup', 'num_duplicate', 'x1']
df_duplicate.head()
df_duplicate.columns
df_duplicate['smpl_id'] = df_duplicate['f'].apply(lambda x: "-".join(x.split("-")[:-1]))
df_duplicate.head()
df_duplicate.tail()
df_duplicate.drop(['f', 'x1'], axis=1, inplace=True)
df_duplicate = df_duplicate.merge(df_descr, 'left', on='smpl_id')
df_duplicate = df_duplicate.dropna()
df_duplicate['subgroup'][df_duplicate.subgroup.isin(['final'])] = 'total'
df_duplicate.protocol.value_counts()

#merge mapped and duplicates
df_FlSt = df_FlSt.merge(df_mapped[['smpl_id', 'subgroup', 'num_mapped']], 'left', on=['smpl_id','subgroup'])
df_FlSt = df_FlSt.merge(df_duplicate[['smpl_id', 'subgroup', 'num_duplicate']], 'left', on=['smpl_id','subgroup'])
df_FlSt.head()
df_FlSt.shape
df_FlSt.dtypes
df_FlSt['perc_mapped'] = df_FlSt['num_mapped'].astype(float)/df_FlSt['num_alignments']
df_FlSt['perc_duplicate'] = df_FlSt['num_duplicate'].astype(float)/df_FlSt['num_alignments']

# a few summaries
def mysumm(x):
    num_alignments = sum(x['num_alignments'])
    perc_mapped = round(sum(x['perc_mapped'])/len(x['perc_mapped']),4)
    perc_duplicate = round(sum(x['perc_duplicate'])/len(x['perc_duplicate']),4)
    return pd.Series ([num_alignments, perc_mapped, perc_duplicate], index = ['num_alignments', 'perc_mapped',
        'perc_duplicate']) 

summ_alingnemt = df_FlSt.groupby(['subgroup', 'protocol']).apply(mysumm)
summ_alingnemt.to_excel("/mnt/xfs1/home/asalomatov/summ_alignment.xls")

### summarize exome coverage


import os
df_descr.head()
df_ex_list = []
for smpl in df_descr['smpl_id']:
    fam = df_descr['fam_id'][df_descr.smpl_id.isin([smpl])].iloc[0]
    rel = df_descr['relationship'][df_descr.smpl_id.isin([smpl])].iloc[0]
    desc = '_'.join([fam, rel])
    print desc
    fname='/mnt/xfs1/scratch/asalomatov/BloodVsSaliva/output/'+smpl+'-D1.final-genes-genes.bed'
    if not os.path.isfile(fname):
        print fname
        continue
    df_exome = pd.read_csv(fname, sep="\t", header=None) 
    df_exome = df_exome[df_exome[0] == 'all']
    print df_exome.head()
    #df_exome.drop([0,1,2, 4, 6, 7], axis=1, inplace=True)
    #df_exome.columns = ['gene_name', 'num_alignments', 'perc']
    df_exome.drop([0,1,3, 5,6,7,8], axis=1, inplace=True)
    df_exome.columns = ['num_alignments', 'perc']
    print df_exome.head()
    df_exome[desc] = 1 - df_exome['perc'].cumsum()
    df_exome = df_exome[df_exome['num_alignments'] <= 1000]
    df_ex_list.append(df_exome)

for f in df_ex_list[1:]:
    print f.shape
    f.drop('perc', axis=1, inplace=True)
    df_exome = df_exome.merge(f, 'left', on='num_alignments')
    print df_exome.shape

df_exome.head()
df_exome = df_exome[df_exome['num_alignments'] <= 400]
sns.set(font_scale=1.8)
df_exome.plot(x='num_alignments',  kind='line')
myplot.set_title('Cumulative Exome Coverage')
myplot.set_ylabel('% of Exome')
myplot.set_xlabel('Depth of Coverage')
plt.savefig("/mnt/xfs1/home/asalomatov/Blood_vs_Saliva_exome_depth.png")
plt.close()

## try a scatter plot
cvrg_20 = df_exome[df_exome['num_alignments'] == 20]
cvrg_20.drop('num_alignments', axis=1, inplace=20)
cvrg_20 = cvrg_20.transpose()
cvrg_20['id'] = cvrg_20.index
cvrg_20.head()
cvrg_20.columns = ['cvrg_perc', 'id']
cvrg_20.reset_index(drop=True, inplace=True)
cvrg_20['protocol'] = cvrg_20.id.apply(lambda x: x.split("_")[0].split("-")[1])
cvrg_20= cvrg_20[~cvrg_20.id.isin(['3-blood_mother'])]
cvrg_20= cvrg_20[~cvrg_20.id.isin(['4-blood_sibling'])]

cvrg_20_melt = pd.melt(cvrg_20, id_vars=['protocol'], value_vars=['cvrg_perc'])
x = cvrg_20_melt['value'][cvrg_20_melt.protocol.isin(['blood']) ]
y = cvrg_20_melt['value'][cvrg_20_melt.protocol.isin(['saliva'])]
plt.scatter(x, y)
plt.xlabel('blood')
plt.ylabel('saliva')
plt.title('Blood vs. Saliva at 20')

sns.set(font_scale=1.8)
cvrg_20.head()
#myplot = sns.barplot(x="id", y="cvrg_perc", data=cvrg_20)
myplot = sns.barplot(x="id", y="cvrg_perc", hue="protocol", data=cvrg_20)
myplot.set_xticklabels(cvrg_20['id'], rotation=90)
myplot.set_title('Mapped')


#myplot = counts_sum.plot(kind='bar', legend=None, title="Number of Alignments")
#sns.barplot(x='protocol', y='num_alignments',  data=df_FlSt)
df_FlSt.head()
df_FlSt = df_FlSt[df_FlSt.subgroup.isin(['total', 'exome'])]
sns.set(font_scale=1.8)
myplot = sns.barplot(x="fam", y="num_alignments", hue="protocol", data=df_FlSt)
myplot.set_title('Number of Alignments')

sns.set(font_scale=1.8)
myplot = sns.barplot(x="fam", y="perc_mapped", hue="protocol", data=df_FlSt)
myplot.set_title('Mapped')

sns.set(font_scale=1.8)
myplot = sns.barplot(x="subgroup", y="perc_duplicate", hue="protocol", data=df_FlSt)
myplot.set_title('Duplicates')
df_FlSt.head()


### now by genes
df_genes = pd.read_csv('/mnt/scratch/asalomatov/data/b37/ensembl-genes-coding_all.bed', sep="\t", header=None) 
df_genes.head()
df_genes.columns = ['chr', 'start', 'end', 'name', 'id']

import os
df_descr.head()
genes_20_blood = set()
genes_20_saliva = set()
for smpl in [df_descr['smpl_id'][0]]:
   # fam = df_descr['fam_id'][df_descr.smpl_id.isin([smpl])].iloc[0]
   # rel = df_descr['relationship'][df_descr.smpl_id.isin([smpl])].iloc[0]
    prot = df_descr['protocol'][df_descr.smpl_id.isin([smpl])].iloc[0]
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
        print cont, sta
        g_list = data.gene_names_at_locus(contig=cont, position=int(sta))
        for g in g_list:
            mygenes.add(g)
    if prot == 'blood':
        genes_20_blood.union(mygenes)
    else:
        genes_20_saliva.union(mygenes)

len(mygenes)
print mygenes
print len(genes_20_blood.intersection(genes_20_saliva))
print len(genes_20_blood.union(genes_20_saliva))
print list(genes_20_blood) + list(genes_20_saliva)
from matplotlib_venn import venn2
venn2([genes_20_blood, genes_20_saliva], set_labels=('genes_20_blood', 'genes_20_saliva') )
plt.savefig("/mnt/xfs1/home/asalomatov/Blood_vs_Saliva_genes_lt_20.png")
plt.close()


### annotate with genes
from pyensembl import EnsemblRelease
data = EnsemblRelease(75)
data.gene_names_at_locus(contig=1, position=100000)
df_exome['gene_name'] = df_exome.apply(lambda x: data.gene_names_at_locus(contig=x['chr'], position=x['start']))
x = df_exome.apply(lambda x: data.gene_names_at_locus(contig=x['chr'], position=(x['start'] + x['end'])/2),
        axis=1)
x

df_exome.head()
df_exome.tail()
df_exome.shape
df_exome.bin.value_counts()
df_exome.isnull().sum()

my_plot.set_xlabel("Customers")
my_plot.set_ylabel("Sales ($)")

x = np.random.normal(size=1000)
y = np.random.normal(size=1000)
plt.scatter(x, y)
plt.plot(x)
sns.distplot(x)
plt.clf()


#plot univariate distributions
x = np.random.normal(size=1000)
print x
sns.distplot(x)

### investigate fam 3-4 blood vs saliva differences
import pandas as pd
import numpy as np
import os

## for NYGC read vcfs to pandas df and summarize genotypes
## for Baylor split big vcf file into families and do same as above
all_list = []
for f in ['1-blood', '2-blood', '3-blood', '4-blood', '1-saliva', '2-saliva', '3-saliva', '4-saliva']:
    fname = "/mnt/scratch/asalomatov/BloodVsSaliva/%s-target.txt" % f
    print fname
    x = pd.read_csv(fname, sep="\t", header=None)
    x.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filt', 'info', 'fmt', 'gt_p1', 'gt_fa', 'gt_mo']
    x.head()
    x['gt_p1'] = x.gt_p1.apply(lambda x: x[:3])
    x['gt_fa'] = x.gt_fa.apply(lambda x: x[:3])
    x['gt_mo'] = x.gt_mo.apply(lambda x: x[:3])
    x['descr'] = f
    x['fam'] = f.split("-")[0]
    x['protocol'] = f.split("-")[1]
    all_list.append(x)

len(all_list)
for i in all_list:
    print i.shape

nyg = pd.concat(all_list)
def varType(row):
    if len(row['ref']) == len(row['alt']):
        return 'snp'
    else:
        return 'indel'
nyg['vartype'] = nyg.apply(varType, axis=1)
nyg.shape
nyg.head()
nyg['gt'] = nyg['gt_p1'] +  '_' + nyg['gt_fa'] + '_' + nyg['gt_mo']
nyg['pos_var'] = nyg['chr'].map(str) + '_' + nyg.pos.map(str)
nyg['pos_allele_var'] = nyg['chr'].map(str) + '_' + nyg.pos.map(str) + '_' +  nyg.ref.map(str)  + '_' +  nyg.alt.map(str) 
nyg['pos_allele_gt_var'] = nyg['chr'].map(str) + '_' + nyg.pos.map(str) + '_' +  nyg.ref.map(str)  + '_' +  nyg.alt.map(str) + '_' +  nyg['gt'] 
nyg['source'] = 'nyg'
nyg = nyg[~nyg.apply(lambda x: "." in x['gt'], axis=1)]

#from functools import partial

def varsumm(x):
    set_bl_pos = set(x['pos_var'][x['protocol']=='blood'])
    set_sl_pos = set(x['pos_var'][x['protocol']=='saliva'])
    set_bl_al = set(x['pos_allele_var'][x['protocol']=='blood'])
    set_sl_al = set(x['pos_allele_var'][x['protocol']=='saliva'])
    set_bl_al_gt = set(x['pos_allele_gt_var'][x['protocol']=='blood'])
    set_sl_al_gt = set(x['pos_allele_gt_var'][x['protocol']=='saliva'])
    var_pos = float(len(set_bl_pos.intersection(set_sl_pos)))/len(set_bl_pos.union(set_sl_pos))
    var_pos_al = float(len(set_bl_al.intersection(set_sl_al)))/len(set_bl_al.union(set_sl_al))
    var_pos_al_gt = float(len(set_bl_al_gt.intersection(set_sl_al_gt)))/len(set_bl_al_gt.union(set_sl_al_gt))
    var_pos_N = len(set_bl_pos.union(set_sl_pos))
    var_pos_al_N = len(set_bl_al.union(set_sl_al))
    var_pos_al_gt_N = len(set_bl_al_gt.union(set_sl_al_gt))
    res = [round(x,2) for x in [var_pos, var_pos_al, var_pos_al_gt]]
    res = res + [var_pos_N]
    return pd.Series (res, index = ['var_pos', 'var_pos_al', 'var_pos_al_gt', 'var_pos_N']) 

nyg.groupby(['fam']).apply(varsumm)
nyg.groupby(['vartype','fam']).apply(varsumm)
nyg.groupby(['fam', 'gt_p1']).apply(varsumm)
x = nyg.groupby(['fam', 'gt_p1']).apply(varsumm)
x[x['var_pos_N'] > 1000]


## Baylor
all_list = []
for f in ['1-blood', '2-blood', '3-blood', '4-blood', '1-saliva', '2-saliva', '3-saliva', '4-saliva']:
    fname = "/mnt/scratch/asalomatov/baylor_qc/%s-baylor.txt" % f
    #fname = "/mnt/scratch/asalomatov/baylor_qc/%s-baylor-indels.txt" % f
    fam = f.split("-")[0]
    prot = f.split("-")[1] 
    if prot == 'saliva' and fam == '3':
        continue
    print fname
    x = pd.read_csv(fname, sep="\t", header=None)
    x.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filt', 'info', 'fmt', 'gt_p1', 'gt_fa', 'gt_mo']
    x.head()
    x['gt_p1'] = x.gt_p1.apply(lambda x: x[:3])
    x['gt_fa'] = x.gt_fa.apply(lambda x: x[:3])
    x['gt_mo'] = x.gt_mo.apply(lambda x: x[:3])
    x['descr'] = f
    x['fam'] = fam
    x['protocol'] = prot
    all_list.append(x)

len(all_list)
for i in all_list:
    print i.shape

for i in nyg.columns:
    if not i in baylor.columns:
        print i

baylor = pd.concat(all_list)
baylor.shape
baylor.head()
baylor['gt'] = baylor['gt_p1'] +  '_' + baylor['gt_fa'] + '_' + baylor['gt_mo']
baylor['pos_var'] = baylor['chr'].map(str) + '_' + baylor.pos.map(str)
baylor['pos_allele_var'] = baylor['chr'].map(str) + '_' + baylor.pos.map(str) + '_' +  baylor.ref.map(str)  + '_' +  baylor.alt.map(str) 
baylor['pos_allele_gt_var'] = baylor['chr'].map(str) + '_' + baylor.pos.map(str) + '_' +  baylor.ref.map(str)  + '_' +  baylor.alt.map(str) + '_' +  baylor['gt'] 
baylor['source'] = 'baylor'
baylor['gt'].value_counts()
baylor = baylor[~baylor['gt'].isin(['0/0_0/0_0/0'])]
baylor['vartype'] = baylor.apply(varType, axis=1)
baylor = baylor[~baylor.apply(lambda x: "." in x['gt'], axis=1)]

## my rerun
all_list = []
for f in ['1', '2', '4']:
    for prot in ['blood', 'saliva']:
        for caller in ['-HC', '-FB', '-PL']:
            #fname = "/mnt/scratch/asalomatov/baylor_qc/rerun/"+prot+"/targ/"+f+caller+"-targ.txt"
            fname = "/mnt/scratch/asalomatov/BloodVsSaliva/rerun/"+prot+"/targ/"+f+caller+"-targ.txt"
            print fname
            x = pd.read_csv(fname, sep="\t", header=None)
            x.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filt', 'info', 'fmt', 'gt_mo', 'gt_fa', 'gt_p1']
            if caller in ['-HC', '-PL']:
                x.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filt', 'info', 'fmt', 'gt_fa', 'gt_mo', 'gt_p1']
            x.head()
            x['gt_p1'] = x.gt_p1.apply(lambda x: x[:3])
            x['gt_fa'] = x.gt_fa.apply(lambda x: x[:3])
            x['gt_mo'] = x.gt_mo.apply(lambda x: x[:3])
            x['descr'] = f
            x['fam'] = f
            x['protocol'] = prot
            all_list.append(x)

len(all_list)
for i in all_list:
    print i.shape

sf_baylor_bam = pd.concat(all_list)
sf_baylor_bam.shape
sf_baylor_bam.columns
sf_baylor_bam.head()
sf_baylor_bam['gt'] = sf_baylor_bam['gt_p1'] +  '_' + sf_baylor_bam['gt_fa'] + '_' + sf_baylor_bam['gt_mo']
sf_baylor_bam['pos_var'] = sf_baylor_bam['chr'].map(str) + '_' + sf_baylor_bam.pos.map(str)
sf_baylor_bam['pos_allele_var'] = sf_baylor_bam['chr'].map(str) + '_' + sf_baylor_bam.pos.map(str) + '_' +  sf_baylor_bam.ref.map(str)  + '_' +  sf_baylor_bam.alt.map(str) 
sf_baylor_bam['pos_allele_gt_var'] = sf_baylor_bam['chr'].map(str) + '_' + sf_baylor_bam.pos.map(str) + '_' +  sf_baylor_bam.ref.map(str)  + '_' +  sf_baylor_bam.alt.map(str) + '_' +  sf_baylor_bam['gt'] 
sf_baylor_bam['source'] = 'sf_baylor_bam'
sf_baylor_bam['gt'].value_counts()
sf_baylor_bam = sf_baylor_bam[~sf_baylor_bam['gt'].isin(['0/0_0/0_0/0'])]
sf_baylor_bam['vartype'] = sf_baylor_bam.apply(varType, axis=1)
sf_baylor_bam = sf_baylor_bam[~sf_baylor_bam.apply(lambda x: "." in x['gt'], axis=1)]
sf_baylor_bam['gt_p1'] = sf_baylor_bam['gt_p1'].apply(lambda x: "/".join([x[0], x[2]]))
sf_baylor_bam['gt_fa'] = sf_baylor_bam['gt_fa'].apply(lambda x: "/".join([x[0], x[2]]))
sf_baylor_bam['gt_mo'] = sf_baylor_bam['gt_mo'].apply(lambda x: "/".join([x[0], x[2]]))
sf_baylor_bam['gt'] = sf_baylor_bam['gt_p1'] +  '_' + sf_baylor_bam['gt_fa'] + '_' + sf_baylor_bam['gt_mo']
sf_baylor_bam['pos_var'] = sf_baylor_bam['chr'].map(str) + '_' + sf_baylor_bam.pos.map(str)
sf_baylor_bam['pos_allele_var'] = sf_baylor_bam['chr'].map(str) + '_' + sf_baylor_bam.pos.map(str) + '_' +  sf_baylor_bam.ref.map(str)  + '_' +  sf_baylor_bam.alt.map(str) 
sf_baylor_bam['pos_allele_gt_var'] = sf_baylor_bam['chr'].map(str) + '_' + sf_baylor_bam.pos.map(str) + '_' +  sf_baylor_bam.ref.map(str)  + '_' +  sf_baylor_bam.alt.map(str) + '_' +  sf_baylor_bam['gt'] 
sf_baylor_bam['source'] = 'sf_baylor_bam'
sf_baylor_bam['gt'].value_counts()
sf_baylor_bam = sf_baylor_bam[~sf_baylor_bam['gt'].isin(['0/0_0/0_0/0'])]

#df['dedup'] = df['pos_allele_gt_var'] + "_" + df['descr']
#df['dedup'].head()
#df = df.drop_duplicates(subset='dedup', take_last=True)

sf_nyg_bam = pd.concat(all_list)
sf_nyg_bam.shape
sf_nyg_bam.head()
sf_nyg_bam['gt'] = sf_nyg_bam['gt_p1'] +  '_' + sf_nyg_bam['gt_fa'] + '_' + sf_nyg_bam['gt_mo']
sf_nyg_bam['pos_var'] = sf_nyg_bam['chr'].map(str) + '_' + sf_nyg_bam.pos.map(str)
sf_nyg_bam['pos_allele_var'] = sf_nyg_bam['chr'].map(str) + '_' + sf_nyg_bam.pos.map(str) + '_' +  sf_nyg_bam.ref.map(str)  + '_' +  sf_nyg_bam.alt.map(str) 
sf_nyg_bam['pos_allele_gt_var'] = sf_nyg_bam['chr'].map(str) + '_' + sf_nyg_bam.pos.map(str) + '_' +  sf_nyg_bam.ref.map(str)  + '_' +  sf_nyg_bam.alt.map(str) + '_' +  sf_nyg_bam['gt'] 
sf_nyg_bam['source'] = 'sf_nyg_bam'
sf_nyg_bam['gt'].value_counts()
sf_nyg_bam = sf_nyg_bam[~sf_nyg_bam['gt'].isin(['0/0_0/0_0/0'])]
sf_nyg_bam['vartype'] = sf_nyg_bam.apply(varType, axis=1)
sf_nyg_bam = sf_nyg_bam[~sf_nyg_bam.apply(lambda x: "." in x['gt'], axis=1)]
sf_nyg_bam['gt_p1'] = sf_nyg_bam['gt_p1'].apply(lambda x: "/".join([x[0], x[2]]))
sf_nyg_bam['gt_fa'] = sf_nyg_bam['gt_fa'].apply(lambda x: "/".join([x[0], x[2]]))
sf_nyg_bam['gt_mo'] = sf_nyg_bam['gt_mo'].apply(lambda x: "/".join([x[0], x[2]]))
sf_nyg_bam['gt'] = sf_nyg_bam['gt_p1'] +  '_' + sf_nyg_bam['gt_fa'] + '_' + sf_nyg_bam['gt_mo']
sf_nyg_bam['pos_var'] = sf_nyg_bam['chr'].map(str) + '_' + sf_nyg_bam.pos.map(str)
sf_nyg_bam['pos_allele_var'] = sf_nyg_bam['chr'].map(str) + '_' + sf_nyg_bam.pos.map(str) + '_' +  sf_nyg_bam.ref.map(str)  + '_' +  sf_nyg_bam.alt.map(str) 
sf_nyg_bam['pos_allele_gt_var'] = sf_nyg_bam['chr'].map(str) + '_' + sf_nyg_bam.pos.map(str) + '_' +  sf_nyg_bam.ref.map(str)  + '_' +  sf_nyg_bam.alt.map(str) + '_' +  sf_nyg_bam['gt'] 
sf_nyg_bam['source'] = 'sf_nyg_bam'
sf_nyg_bam['gt'].value_counts()
sf_nyg_bam = sf_nyg_bam[~sf_nyg_bam['gt'].isin(['0/0_0/0_0/0'])]


#from functools import partial

def varsumm(x):
    set_bl_pos = set(x['pos_var'][x['protocol']=='blood'])
    set_sl_pos = set(x['pos_var'][x['protocol']=='saliva'])
    set_bl_al = set(x['pos_allele_var'][x['protocol']=='blood'])
    set_sl_al = set(x['pos_allele_var'][x['protocol']=='saliva'])
    set_bl_al_gt = set(x['pos_allele_gt_var'][x['protocol']=='blood'])
    set_sl_al_gt = set(x['pos_allele_gt_var'][x['protocol']=='saliva'])
    var_pos = float(len(set_bl_pos.intersection(set_sl_pos)))/len(set_bl_pos.union(set_sl_pos))
    #var_pos_blood = float(len(set_bl_pos.intersection(set_sl_pos)))/len(set_bl_pos)
    #var_pos_saliva = float(len(set_bl_pos.intersection(set_sl_pos)))/len(set_sl_pos)
    var_pos_al = float(len(set_bl_al.intersection(set_sl_al)))/len(set_bl_al.union(set_sl_al))
    #var_pos_al_blood = float(len(set_bl_al.intersection(set_sl_al)))/len(set_bl_al)
    #var_pos_al_saliva = float(len(set_bl_al.intersection(set_sl_al)))/len(set_sl_al)
    var_pos_al_gt = float(len(set_bl_al_gt.intersection(set_sl_al_gt)))/len(set_bl_al_gt.union(set_sl_al_gt))
    #var_pos_al_gt_blood = float(len(set_bl_al_gt.intersection(set_sl_al_gt)))/len(set_bl_al_gt)
    #var_pos_al_gt_saliva = float(len(set_bl_al_gt.intersection(set_sl_al_gt)))/len(set_sl_al_gt)
    var_pos_N = len(set_bl_pos.union(set_sl_pos))
    var_pos_al_N = len(set_bl_al.union(set_sl_al))
    var_pos_al_gt_N = len(set_bl_al_gt.union(set_sl_al_gt))
    #res = [round(x,2) for x in [var_pos,var_pos_blood,var_pos_saliva, var_pos_al,var_pos_al_blood,var_pos_al_saliva,
    #    var_pos_al_gt, var_pos_al_gt_blood, var_pos_al_gt_saliva]]
    res = [round(x,2) for x in [var_pos, var_pos_al, var_pos_al_gt]]
    res = res + [var_pos_N, var_pos_al_N, var_pos_al_gt_N]
    #return pd.Series (res, index = ['var_pos/union','var_pos/blood','var_pos/saliva', 'var_pos_al/union',
    #    'var_pos_al/blood', 'var_pos_al/saliva', 'var_pos_al_gt/union', 'var_pos_al_gt/blood', 'var_pos_al_gt/saliva', 'var_pos_N']) 
    return pd.Series (res, index = ['var_pos/union', 'var_pos_al/union', 'var_pos_al_gt/union', 'var_pos_N',
        'var_pos_al_N', 'var_pos_al_gt_N']) 

def myfunc(x):
    if len(x) == 3:
        return x+'0'
    else:
        return x

def protsumm(x, labels):
    label1 = labels[0]
    label2 = labels[1]
    c1 = x['source'] == label1
    c2 = x['source'] == label2
    set_bl_pos = set(x['pos_var'][c1])
    set_sl_pos = set(x['pos_var'][c2])
    set_bl_al = set(x['pos_allele_var'][c1])
    set_sl_al = set(x['pos_allele_var'][c2])
    set_bl_al_gt = set(x['pos_allele_gt_var'][c1])
    set_sl_al_gt = set(x['pos_allele_gt_var'][c2])
    var_pos = float(len(set_bl_pos.intersection(set_sl_pos)))/len(set_bl_pos.union(set_sl_pos))
    var_pos_bl = float(len(set_bl_pos.intersection(set_sl_pos)))/len(set_bl_pos)
    var_pos_sl = float(len(set_bl_pos.intersection(set_sl_pos)))/len(set_sl_pos)
    var_pos_al = float(len(set_bl_al.intersection(set_sl_al)))/len(set_bl_al.union(set_sl_al))
    var_pos_al_bl = float(len(set_bl_al.intersection(set_sl_al)))/len(set_bl_al)
    var_pos_al_sl = float(len(set_bl_al.intersection(set_sl_al)))/len(set_sl_al)
    var_pos_al_gt = float(len(set_bl_al_gt.intersection(set_sl_al_gt)))/len(set_bl_al_gt.union(set_sl_al_gt))
    var_pos_al_gt_bl = float(len(set_bl_al_gt.intersection(set_sl_al_gt)))/len(set_bl_al_gt)
    var_pos_al_gt_sl = float(len(set_bl_al_gt.intersection(set_sl_al_gt)))/len(set_sl_al_gt)
    var_pos_N = len(set_bl_pos.union(set_sl_pos))
    var_pos_al_N = len(set_bl_al.union(set_sl_al))
    var_pos_al_gt_N = len(set_bl_al_gt.union(set_sl_al_gt))
    res = [str(round(x,2)) for x in [var_pos_bl,var_pos,var_pos_sl, var_pos_al_bl,var_pos_al,var_pos_al_sl,var_pos_al_gt_bl,
        var_pos_al_gt,var_pos_al_gt_sl]]
    res0 = map(myfunc, res)
    res1 = ['/'.join(res0[:3]), '/'.join(res0[3:6]),'/'.join(res0[6:])]
    res1 = res1 + [var_pos_N, var_pos_al_N, var_pos_al_gt_N]
    return pd.Series (res1, index = ['/'.join([label1,'var_pos',label2]), '/'.join([label1,'var_pos_al',label2]),
        '/'.join([label1,'var_pos_al_gt',label2]), 'var_pos_N', 'var_pos_al_N', 'var_pos_al_gt_N']) 


baylor.groupby(['fam']).apply(varsumm)
baylor.groupby(['vartype','fam']).apply(varsumm)
baylor.groupby(['fam', 'gt_p1']).apply(varsumm)
x = baylor.groupby(['fam', 'gt']).apply(varsumm)
x[x['var_pos_N'] > 1000]


baylor.groupby(['vartype','fam']).apply(varsumm)
baylor.groupby(['fam', 'gt_p1']).apply(varsumm)
x = baylor.groupby(['fam', 'gt']).apply(varsumm)
x[x['var_pos_N'] > 1000]

### concat nyg and baylor
df = pd.concat([nyg, baylor])
df.shape
sf_nyg_bam.fam.value_counts()
df_all = pd.concat([df, sf_baylor_bam, sf_nyg_bam])
df_all = df_all[df_all['fam'] != '3']
df_allCall = pd.concat([df, sf_baylor_bam, sf_nyg_bam])
df_allCall = df_allCall[df_allCall['fam'] != '3']
df['dedup'] = df['pos_allele_gt_var'] + "_" + df['descr']
df['dedup'].head()
df = df.drop_duplicates(subset='dedup', take_last=True)

df.groupby(['source','fam']).apply(varsumm)
df_all.groupby(['source','fam']).apply(varsumm)
df_all.groupby(['source','fam']).apply(varsumm).to_excel('/mnt/xfs1/home/asalomatov/nov10_blood_vs_saliva_1.xls')
df_allCall.groupby(['source','fam']).apply(varsumm)
df_allCall.groupby(['source','fam']).apply(varsumm).to_excel('/mnt/xfs1/home/asalomatov/nov10_blood_vs_saliva_1.xls')
df.groupby(['source','fam']).apply(varsumm).to_excel('/mnt/xfs1/home/asalomatov/blood_vs_saliva_1.xls')
df.groupby(['source','vartype','fam']).apply(varsumm)
df.groupby(['source','vartype','fam']).apply(varsumm).to_excel('/mnt/xfs1/home/asalomatov/blood_vs_saliva_2.xls') 


df.groupby(['protocol', 'fam']).apply(protsumm)
df_all.groupby(['protocol', 'fam']).apply(protsumm, ('baylor', 'nyg'))
df_all.groupby(['protocol', 'fam']).apply(protsumm, ('sf_nyg_bam', 'nyg'))
df_all.groupby(['protocol', 'fam']).apply(protsumm, ('sf_nyg_bam', 'sf_baylor_bam'))
df_all.groupby(['protocol', 'fam']).apply(protsumm, ('sf_nyg_bam', 'baylor'))
df_all.groupby(['protocol', 'fam']).apply(protsumm, ('sf_nyg_bam', 'nyg')).to_excel('/mnt/xfs1/home/asalomatov/nov10_nyg_vs_sfHC.xls') 

df_all.groupby(['protocol', 'fam']).apply(protsumm, ('sf_baylor_bam', 'baylor'))
df_all.groupby(['protocol', 'fam']).apply(protsumm, ('sf_baylor_bam', 'baylor')).to_excel('/mnt/xfs1/home/asalomatov/nov10_baylor_vs_sfHC.xls') 
df_all.groupby(['protocol', 'fam', 'vartype']).apply(protsumm, ('sf_baylor_bam', 'baylor'))
df_all.groupby(['protocol', 'fam', 'vartype']).apply(protsumm, ('sf_baylor_bam', 'baylor')).to_excel('/mnt/xfs1/home/asalomatov/nov10_baylor_vs_sfHC_vartype.xls') 

df_allCall.groupby(['protocol', 'fam']).apply(protsumm, ('sf_baylor_bam', 'baylor'))
df_allCall.groupby(['protocol', 'fam']).apply(protsumm, ('sf_baylor_bam','baylor')).to_excel('/mnt/xfs1/home/asalomatov/nov10_baylor_vs_sfHCFBPL.xls') 
df_allCall.groupby(['protocol', 'fam', 'vartype']).apply(protsumm, ('sf_baylor_bam', 'baylor'))
df_allCall.groupby(['protocol', 'fam', 'vartype']).apply(protsumm, ('sf_baylor_bam',
    'baylor')).to_excel('/mnt/xfs1/home/asalomatov/nov10_baylor_vs_sfHCFBPL_vartype.xls') 

df_all.groupby(['protocol', 'fam']).apply(protsumm, ('sf_baylor_bam', 'sf_nyg_bam'))
df_all.groupby(['protocol', 'fam']).apply(protsumm, ('sf_baylor_bam', 'sf_nyg_bam')).to_excel('/mnt/xfs1/home/asalomatov/nov10_sfHCbay_vs_sfHCnyg.xls') 
df_all.groupby(['protocol', 'fam', 'vartype']).apply(protsumm, ('sf_baylor_bam', 'sf_nyg_bam'))
df_all.groupby(['protocol', 'fam', 'vartype']).apply(protsumm, ('sf_baylor_bam', 'sf_nyg_bam')).to_excel('/mnt/xfs1/home/asalomatov/nov10_sfHCbay_vs_sfHCnyg_vartype.xls') 

df_allCall.groupby(['protocol', 'fam']).apply(protsumm, ('sf_baylor_bam', 'sf_nyg_bam'))
df_allCall.groupby(['protocol', 'fam']).apply(protsumm, ('sf_baylor_bam',
    'sf_nyg_bam')).to_excel('/mnt/xfs1/home/asalomatov/nov10_sfHCFBPLbay_vs_sfHCFBPLnyg.xls') 
df_allCall.groupby(['protocol', 'fam', 'vartype']).apply(protsumm, ('sf_baylor_bam', 'sf_nyg_bam'))
df_allCall.groupby(['protocol', 'fam', 'vartype']).apply(protsumm, ('sf_baylor_bam',
    'sf_nyg_bam')).to_excel('/mnt/xfs1/home/asalomatov/nov10_sfHCFBPLbay_vs_sfHCFBPLnyg_vartype.xls') 

df_all.groupby(['protocol', 'fam']).apply(protsumm1)
df.groupby(['protocol', 'fam']).apply(protsumm).to_excel('/mnt/xfs1/home/asalomatov/nyg_vs_baylor_1.xls') 
df.groupby(['protocol', 'fam', 'vartype']).apply(protsumm)
df.groupby(['protocol', 'fam', 'vartype']).apply(protsumm).to_excel('/mnt/xfs1/home/asalomatov/nyg_vs_baylor_2.xls') 

df_all['gt_p1'].value_counts()
df_allCall['gt_p1'].value_counts()
df_all['source'].value_counts()
x = df_all['fam'].value_counts()
x
'sf_baylor_bam' in x
'sf_nyg_bam' in x
'nyg' in x
'baylor' in x

label1 = 'sf_nyg_bam'
label2 = 'sf_baylor_bam'
c1 = df_all['source'] == label1
c2 = df_all['source'] == label2
sum(c1)
sum(c2)
set_bl_pos = set(df_all['pos_var'][c1])
set_sl_pos = set(df_all['pos_var'][c2])
var_pos = float(len(set_bl_pos.intersection(set_sl_pos)))/len(set_bl_pos.union(set_sl_pos))
var_pos

## look at examples of discordant variants
#for nyg

c1 = df['source'] == 'baylor'
c2 = df['vartype'] == 'snp'
c3 = df['fam'] == '1'
c4 = df['protocol'] == 'blood'
sum(c1 & c2 & c3 & c4)
sum(c1 & c2 & c3 & ~c4)

set_bl_pos = set(df['pos_var'][c1 & c2 & c3 & c4])
set_sl_pos = set(df['pos_var'][c1 & c2 & c3 & ~c4])
set_pos_isec = set_bl_pos.intersection(set_sl_pos)
set_bl_al = set(df['pos_allele_var'][c1 & c2 & c3 & c4])
set_sl_al = set(df['pos_allele_var'][c1 & c2 & c3 & ~c4])
set_al_isec = set_bl_al.intersection(set_sl_al)
set_bl_al_gt = set(df['pos_allele_gt_var'][c1 & c2 & c3 & c4])
set_sl_al_gt = set(df['pos_allele_gt_var'][c1 & c2 & c3 & ~c4])
set_al_gt_isec = set_bl_al_gt.intersection(set_sl_al_gt)

mysubset = df[c1 & c2 & c3]
mysubset['qual'].describe()
df['qual'].describe()
df['qual'][c1].describe()
df['qual'][~c1].describe()

mysubset['pos_TF'] = mysubset['pos_var'].apply(lambda a: a in set_pos_isec)
mysubset['pos_al_TF'] = mysubset['pos_allele_var'].apply(lambda a: a in set_al_isec)
mysubset['pos_gt_TF'] = mysubset['pos_allele_gt_var'].apply(lambda a: a in set_al_gt_isec)
mysubset['pos_TF'].value_counts()
mysubset['pos_al_TF'].value_counts()
mysubset['pos_gt_TF'].value_counts()

cond1 = mysubset['pos_TF'] & ~mysubset['pos_al_TF']
mysubset[['pos_allele_var','qual', 'filt', 'info', 'gt', 'fam', 'protocol']][cond1]

cond2 = mysubset['pos_al_TF'] & ~mysubset['pos_gt_TF']
c5 = (mysubset['gt_p1'] == '0/1') & (mysubset['chr'] == '1') 
fam3_mismatch_gt =  mysubset[['pos_allele_gt_var','qual', 'filt', 'info', 'gt', 'fam', 'protocol']][cond2 & c5]
fam3_mismatch_gt.shape
fam3_mismatch_gt.sort(['pos_allele_gt_var']).tail(20)
fam3_mismatch_gt['qual'].mean()

df[(df['pos_allele_var'] == '10_100018844_G_A') & (df['fam'] == '4')]
for i in df['info'][(df['pos_allele_var'] == '10_100018844_G_A') & (df['fam'] == '4')]:
    print i



