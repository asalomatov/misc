
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
    fname='/mnt/xfs1/scratch/asalomatov/BloodVsSaliva/output/'+smpl+'-D1.final-exome-cvrg.txt'
    if not os.path.isfile(fname):
        print fname
        continue
    df_exome = pd.read_csv(fname, sep="\t", header=None) 
    df_exome = df_exome[df_exome[0] == 'all']
    df_exome.drop([0, 2, 3, 5, 6], axis=1, inplace=True)
    df_exome.columns = ['num_alignments', 'perc']
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

