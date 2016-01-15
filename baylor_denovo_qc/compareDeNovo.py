import pandas as pd
import numpy as np
import os

geneDX = pd.DataFrame({'chr':['20', '4', '14', '4', '14'],'pos':['31021623', '501248', '21854332', '501248', '21854332'],'ref':['-', 'G', 'G', 'G', 'G'],'alt':['AAAAGAAGCCCCG', 'A', 'T', 'A', 'T'],'descr':['frameshift', 'missense', 'missense', 'missense', 'missense'],'gene':['ASXL1', 'CHD8',
        'PIGG', 'CHD8', 'PIGG'],'fam':['fam2', 'fam3', 'fam3', 'fam3', 'fam3'],'memb':['p1', 'p1', 'p1', 'p1',
            'p1'],'prot':['?', 'blood', 'blood', 'saliva', 'saliva']}, index=range(5))
geneDX['source'] = 'geneDX'


df_bay_mendelian = pd.read_csv("/mnt/scratch/asalomatov/baylor_qc/baylor_mendelian.txt", header=None, sep="\t")
df_bay_mendelian.head()
df_bay_mendelian.isnull().sum()
df_bay_mendelian.dropna(axis=1, inplace=True)
df_bay_mendelian.columns = ['chr', 'pos', 'ref', 'alt', 'descr', 'gene', 'fam', 'memb', 'prot']
df_bay_mendelian.head()
df_bay_mendelian['source'] = 'baylor_mendelian'

df_bay_codified = pd.read_csv("/mnt/scratch/asalomatov/baylor_qc/baylor_codified.txt", header=None, sep="\t")
df_bay_codified.head()
df_bay_codified.isnull().sum()
df_bay_codified.dropna(axis=1, inplace=True)
df_bay_codified.columns = ['chr', 'pos', 'ref', 'alt', 'descr', 'gene', 'fam', 'memb', 'prot']
df_bay_codified.head()
df_bay_codified['source'] = 'baylor_codified'

baylor_denovo = pd.concat([df_bay_mendelian, df_bay_codified])
baylor_denovo['vartype'] = baylor_denovo.apply(varType, axis=1)
baylor_denovo.head(50)
baylor_denovo.tail()
baylor_denovo.pos_var
#from natsort import natsorted
#natsorted(baylor_denovo.pos_var)

baylor_geneDX_denovo = pd.concat([df_bay_mendelian, df_bay_codified, geneDX])
baylor_geneDX_denovo['fam'] = baylor_geneDX_denovo['fam'].apply(lambda x: x[-1])
baylor_geneDX_denovo = pd.concat([baylor_geneDX_denovo, sf_indel_denovo])
baylor_geneDX_denovo['fam_prot_pos_var'] = baylor_geneDX_denovo['fam'] + '_' + baylor_geneDX_denovo['prot'] + '_' + baylor_geneDX_denovo['chr'].map(str) + '_' + baylor_geneDX_denovo.pos.map(str)
baylor_geneDX_denovo['fam_memb_pos_var'] = baylor_geneDX_denovo['fam'] + '_' + baylor_geneDX_denovo['memb'] + '_' + baylor_geneDX_denovo['chr'].map(str) + '_' + baylor_geneDX_denovo.pos.map(str)
baylor_geneDX_denovo['fam_prot_pos_var_al'] = baylor_geneDX_denovo['fam'] + '_' + baylor_geneDX_denovo['prot'] + '_' + baylor_geneDX_denovo['chr'].map(str) + '_' + baylor_geneDX_denovo.pos.map(str) + '_' + baylor_geneDX_denovo.ref + '_' + baylor_geneDX_denovo.alt
baylor_geneDX_denovo['vartype'] = baylor_geneDX_denovo.apply(varType, axis=1)
baylor_geneDX_denovo.head()
baylor_geneDX_denovo.source.value_counts()
baylor_geneDX_denovo.prot.value_counts()
baylor_geneDX_denovo[baylor_geneDX_denovo.vartype.isin(['snp'])]
baylor_geneDX_denovo.head()
baylor_geneDX_denovo.shape
baylor_geneDX_denovo_uniq = baylor_geneDX_denovo.drop_duplicates(subset='fam_memb_pos_var', take_last=True)
baylor_geneDX_denovo_uniq.shape

blrM_bl=[]
blrM_sa=[]
blrC_bl=[]
blrC_sa=[]
nyg_bl=[]
nyg_sa=[]
sf_bl=[]
sf_sa=[]
sf_bl_indel=[]
sf_sa_indel=[]
gdx_bl=[]
gdx_sa=[]
#for i, row in baylor_geneDX_denovo[baylor_geneDX_denovo.vartype.isin(['snp'])].iterrows():
for i, row in baylor_geneDX_denovo_uniq.iterrows():
    c1 = (baylor_geneDX_denovo['source'] == 'baylor_mendelian')
    c2 = baylor_geneDX_denovo['fam_memb_pos_var'] == row['fam_memb_pos_var']
    c3 = baylor_geneDX_denovo['prot'] == 'blood'
    c4 = baylor_geneDX_denovo['prot'] == 'saliva'
    blrM_bl.append(sum(c1 & c2 & c3))
    blrM_sa.append(sum(c1 & c2 & c4))
    c1 = (baylor_geneDX_denovo['source'] == 'baylor_codified')
    c2 = baylor_geneDX_denovo['fam_memb_pos_var'] == row['fam_memb_pos_var']
    c3 = baylor_geneDX_denovo['prot'] == 'blood'
    c4 = baylor_geneDX_denovo['prot'] == 'saliva'
    blrC_bl.append(sum(c1 & c2 & c3))
    blrC_sa.append(sum(c1 & c2 & c4))
    c1 = (baylor_geneDX_denovo['source'] == 'geneDX')
    c2 = baylor_geneDX_denovo['fam_memb_pos_var'] == row['fam_memb_pos_var']
    c3 = baylor_geneDX_denovo['prot'] == 'blood'
    c4 = baylor_geneDX_denovo['prot'] == 'saliva'
    gdx_bl.append(sum(c1 & c2 & c3))
    gdx_sa.append(sum(c1 & c2 & c4))
    c1 = (baylor_geneDX_denovo['source'] == 'sf')
    c2 = baylor_geneDX_denovo['fam_memb_pos_var'] == row['fam_memb_pos_var']
    c3 = baylor_geneDX_denovo['prot'] == 'blood'
    c4 = baylor_geneDX_denovo['prot'] == 'saliva'
    sf_bl.append(sum(c1 & c2 & c3))
    sf_sa.append(sum(c1 & c2 & c4))
    c1 = df_nyg['prot'] == 'blood'
    c2 = df_nyg['prot'] == 'saliva'
    nyg_bl.append(sum(row['fam_prot_pos_var'] == df_nyg['fam_prot_pos_var']))
    nyg_sa.append(sum(row['fam_prot_pos_var'] == df_nyg['fam_prot_pos_var']))

baylor_geneDX_denovo.shape
len(center_prot)
baylor_geneDX_denovo_uniq['BlrM-bl'] = blrM_bl
baylor_geneDX_denovo_uniq['BlrM-sa'] = blrM_sa
baylor_geneDX_denovo_uniq['BlrC-bl'] = blrC_bl
baylor_geneDX_denovo_uniq['BlrC-sa'] = blrC_sa
baylor_geneDX_denovo_uniq['nyg-bl'] = nyg_bl
baylor_geneDX_denovo_uniq['nyg-sa'] = nyg_sa
baylor_geneDX_denovo_uniq['sf-bl'] = sf_bl
baylor_geneDX_denovo_uniq['sf-sa'] = sf_sa
#baylor_geneDX_denovo_uniq['sf-bl-ind'] = sf_bl_indel
#baylor_geneDX_denovo_uniq['sf-sa-ind'] = sf_sa_indel
baylor_geneDX_denovo_uniq['gdx-bl'] = gdx_bl
baylor_geneDX_denovo_uniq['gdx-sa'] = gdx_sa

baylor_geneDX_denovo_uniq[['chr', 'pos', 'ref', 'alt', 'descr', 'gene', 'fam', 'memb', 'BlrM-bl', 'BlrM-sa','BlrC-bl',
    'BlrC-sa','gdx-bl', 'gdx-sa','sf-bl', 'sf-sa','nyg-bl','nyg-sa']][baylor_geneDX_denovo_uniq.vartype.isin(['snp'])]
baylor_geneDX_denovo_uniq[['chr', 'pos', 'ref', 'alt', 'descr', 'gene', 'fam', 'memb', 'BlrM-bl', 'BlrM-sa','BlrC-bl', 'BlrC-sa','gdx-bl', 'gdx-sa','sf-bl',
    'sf-sa','nyg-bl','nyg-sa']][baylor_geneDX_denovo_uniq.vartype.isin(['snp'])].to_excel('/mnt/xfs1/home/asalomatov/nov24_denovo_snp.xls')

baylor_geneDX_denovo_uniq[['chr', 'pos', 'ref', 'alt', 'descr', 'gene', 'fam', 'memb', 'BlrM-bl', 'BlrM-sa','BlrC-bl',
    'BlrC-sa','gdx-bl', 'gdx-sa','sf-bl', 'sf-sa','nyg-bl','nyg-sa']][~baylor_geneDX_denovo_uniq.vartype.isin(['snp'])]
baylor_geneDX_denovo_uniq[['chr', 'pos', 'ref', 'alt', 'descr', 'gene', 'fam', 'memb', 'BlrM-bl', 'BlrM-sa','BlrC-bl',
    'BlrC-sa','gdx-bl', 'gdx-sa','sf-bl', 'sf-sa','nyg-bl','nyg-sa']][~baylor_geneDX_denovo_uniq.vartype.isin(['snp'])].to_excel('/mnt/xfs1/home/asalomatov/nov24_denovo_indel.xls')

df_n
###DNM Filter applied to Baylor calls
df_list = []
for f in [1,2,3,4]:
    for p in ['bl', 'sa']:
        if f==3 and p == 'sa': continue
        file_name = str(f)+'-'+p+'-dnm.csv'
        df = pd.read_csv("/mnt/scratch/asalomatov/baylor_qc/dnm_filter/output/"+file_name, header=False,
            sep=",")
        df.columns = ['desc', 'chr', 'pos', 'score']
        print df.head()
        df_list.append(df)

len(df_list)
for i in df_list:
    print i.shape

df_bay = pd.concat(df_list)
df_bay.head()
df_bay['fam'] = df_bay['desc'].apply(lambda x: x.split('-')[0])
df_bay['prot'] = df_bay['desc'].apply(lambda x: x.split('-')[1])
df_bay.head()
df_bay.chr.value_counts()
df_bay.isnull().sum()
df_bay['pos_var'] = df_bay['chr'].map(str) + '_' + df_bay.pos.map(str)
df_bay['source'] = 'baylor'

df_list = []
for f in [1,2,3,4]:
    for p in ['bl', 'sa']:
        if f==3 and p == 'sa': continue
        file_name = str(f)+'-'+p+'-dnm.csv'
        df = pd.read_csv("/mnt/scratch/asalomatov/baylor_qc/rerun/dnm_filter/output/"+file_name, header=None,
            sep=",")
        df.columns = ['desc', 'chr', 'pos', 'score']
        print df.head()
        df_list.append(df)

len(df_list)
for i in df_list:
    print i.shape

df_rerun = pd.concat(df_list)
df_rerun.head()
df_rerun['fam'] = df_rerun['desc'].apply(lambda x: x.split('-')[0])
df_rerun['prot'] = df_rerun['desc'].apply(lambda x: x.split('-')[1])
df_rerun.head()
df_rerun.chr.value_counts()
df_rerun.isnull().sum()
df_rerun['pos_var'] = df_rerun['chr'].map(str) + '_' + df_rerun.pos.map(str)
df_rerun['source'] = 'rerun'


def myfunc(x, char='0'):
    if len(x) == 3:
        return x+char
    else:
        return x

def protsumm(x, extra_args):
    label1 = extra_args[0]
    label2 = extra_args[1]
    field = extra_args[2]
    score_cutoff = extra_args[3]
    c1 = (x[field] == label1) & (x['score'] > score_cutoff)
    c2 = (x[field] == label2)  & (x['score'] > score_cutoff)
    set_lab1 = set(x['pos_var'][c1])
    set_lab2 = set(x['pos_var'][c2])
    var_pos = float(len(set_lab1.intersection(set_lab2)))/len(set_lab1.union(set_lab2))
    var_pos_lab1 = float(len(set_lab1.intersection(set_lab2)))/len(set_lab1)
    var_pos_lab2 = float(len(set_lab1.intersection(set_lab2)))/len(set_lab2)
    var_pos_N = len(set_lab1.union(set_lab2))
    res = [str(round(x,2)) for x in [var_pos_lab1,var_pos,var_pos_lab2]]
    res0 = map(myfunc, res)
    res1 = ['/'.join(res0[:3])]
    #res1 = ['/'.join(res0[:3]), '/'.join(res0[3:6]),'/'.join(res0[6:])]
    res1 = res1 + [var_pos_N]
    return pd.Series (res1, index = ['/'.join([label1,'var_pos',label2]), 'var_pos_N']) 

def protsummNum(x, extra_args):
    label1 = extra_args[0]
    label2 = extra_args[1]
    field = extra_args[2]
    score_cutoff = extra_args[3]
    c1 = (x[field] == label1) & (x['score'] > score_cutoff)
    c2 = (x[field] == label2)  & (x['score'] > score_cutoff)
    set_lab1 = set(x['pos_var'][c1])
    set_lab2 = set(x['pos_var'][c2])
    var_pos = len(set_lab1.intersection(set_lab2))
    var_pos_lab1 = len(set_lab1)
    var_pos_lab2 = len(set_lab2)
    var_pos_N = len(set_lab1.union(set_lab2))
    #res = [str(round(x,2)) for x in [var_pos_lab1,var_pos,var_pos_lab2]]
    res = map(str, [var_pos_lab1,var_pos,var_pos_lab2])
    #res0 = [myfunc(i, ' ') for i in res]
    res0 = res
    res1 = ['/'.join(res0[:3])]
    #res1 = ['/'.join(res0[:3]), '/'.join(res0[3:6]),'/'.join(res0[6:])]
    res1 = res1 + [var_pos_N]
    return pd.Series (res1, index = ['/'.join([label1,'var_pos',label2]), 'var_pos_N']) 

df_rerun.groupby(['fam']).apply(protsumm, ('bl', 'sa', 'prot',  .5))

rerun_bl_vs_sa_list = []
for sc in [0, .5, .6, .7, .8, .9, .95, .97]:
    print sc
    mysumm =  df_rerun.groupby(['fam']).apply(protsumm, ('bl', 'sa', 'prot',  sc))
    mysumm['score'] = sc
    rerun_bl_vs_sa_list.append(mysumm)

rerun_bl_vs_sa_denovo = pd.concat(rerun_bl_vs_sa_list)
rerun_bl_vs_sa_denovo
rerun_bl_vs_sa_denovo.to_excel('/mnt/xfs1/home/asalomatov/nov16_rerun_bl_vs_sa.xls')
df_rerun.groupby(['fam']).apply(protsumm, ('bl', 'sa')).to_excel('/mnt/xfs1/home/asalomatov/nov16_ueeeeeeeeeeeee.xls') 



baylor_bl_vs_sa_list = []
for sc in [0, .5, .6, .7, .8, .9, .95, .97]:
    print sc
    mysumm =  df_bay.groupby(['fam']).apply(protsumm, ('bl', 'sa', 'prot',  sc))
    mysumm['score'] = sc
    baylor_bl_vs_sa_list.append(mysumm)

baylor_bl_vs_sa_denovo = pd.concat(baylor_bl_vs_sa_list)
baylor_bl_vs_sa_denovo
baylor_bl_vs_sa_denovo.to_excel('/mnt/xfs1/home/asalomatov/nov16_baylor_bl_vs_sa.xls')

df_all = pd.concat([df_rerun, df_bay])
df_all.groupby(['source', 'fam']).apply(protsumm, ('bl', 'sa', 'prot',  .9))
df_all.groupby(['prot', 'fam']).apply(protsumm, ('rerun', 'baylor', 'source',  .95))

temp_list = []
for sc in [0, .5, .6, .7, .8, .9, .95, .97]:
    print sc
    mysumm =  df_all.groupby(['prot', 'fam']).apply(protsumm, ('rerun', 'baylor', 'source',  sc))
    mysumm['score'] = sc
    temp_list.append(mysumm)
len(temp_list)


rerun_vs_baylor_denovo = pd.concat(temp_list)
rerun_vs_baylor_denovo
rerun_vs_baylor_denovo.to_excel('/mnt/xfs1/home/asalomatov/nov16_rerun_vs_baylor_denovo.xls')




### use variants to parse vcf files from our rerun on baylor's bams

import sys
sys.path.insert(0, '/mnt/xfs1/home/asalomatov/projects/variants/variants')
import variants
import ped
import func
import features_vcf as fv
import pandas as pd
import features_vcf
import numpy as np
import os

myvars = variants.Variants('/mnt/scratch/asalomatov/baylor_qc/simons_pvcf_snps.annotated.onTarget.vcf', '1')

myvars.variants.vartype.value_counts()
myvars.variants.columns
df_temp['score'][df_temp.ID.isnull()].mean()
df_temp['score'][~df_temp.ID.isnull()].mean()

run readVarsDnm.py

len(temp_list_snp)
len(temp_list_indel)
for i in temp_list_snp:
    print i.shape

for i in temp_list_indel:
    print i.shape

#temp_list_snp_reun = temp_list_snp
#temp_list_indel_reun = temp_list_indel
df_re_snp = pd.concat(temp_list_snp)
df_re_indel = pd.concat(temp_list_indel)
df_re_snp.shape
df_re_indel.shape
df_re_indel.vartype.value_counts()
df_re_indel.family_id.value_counts()
df_re_indel.columns
df_re_snp.columns
print list(df_re_snp.columns)
df_re_snp['fam_prot_pos_var'] = df_re_snp['fam'].map(str) + '_' + df_re_snp['prot'] + '_' + df_re_snp['CHROM'].map(str) + '_' + df_re_snp.POS.map(str)
df_re_indel['fam_prot_pos_var'] = df_re_indel['fam'].map(str) + '_' + df_re_indel['prot'] + '_' + df_re_indel['CHROM'].map(str) + '_' + df_re_indel.POS.map(str)
sum((df_re_snp['p1_gt_type'] >= 1) & (df_re_snp['fa_gt_type'] <= 1) & (df_re_snp['mo_gt_type'] <= 1))
sum((df_re_snp['p1_gt_type'] >= 1) & (df_re_snp['fa_gt_type'] <= 0) & (df_re_snp['mo_gt_type'] <= 0))

df_re_snp['info_dbNSFP_1000Gp3_AF'].isnull().sum()
df_re_snp['info_dbNSFP_1000Gp3_AF'].describe()
print list(df_re_snp.columns)
df_re_snp.groupby(['fam']).apply(protsumm, ('bl', 'sa', 'prot',  .75))
df_re_snp.groupby(['fam']).apply(protsummNum, ('bl', 'sa', 'prot',  .75))
df_re_snp.groupby(['prot', 'fam']).apply(protsumm, ('rerun', 'baylor', 'source',  .95))

### SNPs with dnmscore > .75 and rs IDs
df_re_snp.columns
df_re_snp_rs_75 = df_re_snp[['family_id', 'CHROM', 'POS', 'ID',
    'REF','ALT','QUAL','format_p1_GT','format_fa_GT','format_mo_GT','score']][(df_re_snp['score']>=.5) & (df_re_snp.ID.isnull())]
df_re_snp_rs_75.shape
df_re_snp_rs_75
df_re_snp_rs_75.to_excel('/mnt/xfs1/home/asalomatov/nov18_noRS_dnm50.xls')

my_list = []
for sc in [0,.3, .4, .5, .6, .7, .75]:
    print sc
    mysumm =  df_re_snp.groupby(['fam']).apply(protsummNum, ('blood', 'saliva', 'prot',  sc))
    mysumm['score'] = sc
    my_list.append(mysumm)

len(my_list)


rerun_vs_baylor_denovo = pd.concat(my_list) rerun_vs_baylor_denovo rerun_vs_baylor_denovo.to_excel('/mnt/xfs1/home/asalomatov/nov16_rerun_vs_baylor_denovo.xls') ### if desired remove snps with rs IDs df_re_snp =  df_re_snp[df_re_snp.ID.isnull()] df_re_snp.shape ###apply hard filter(Ash) df_re_snp['p1_gt_type'].value_counts() c_GT = (df_re_snp['p1_gt_type'] >= 1) & (df_re_snp['fa_gt_type'] <= 0) & (df_re_snp['mo_gt_type'] <= 0) c_AF = (df_re_snp['info_dbNSFP_1000Gp3_AF'] < .001) | (df_re_snp['info_dbNSFP_ExAC_AF'] < .001) | (df_re_snp['info_dbNSFP_ExAC_AF'].isnull()) c_FILTER = (df_re_snp['FILTER'] == 'PASS') c_DP = (df_re_snp['format_p1_DP'] >= 8) & (df_re_snp['format_fa_DP'] >= 8) & (df_re_snp['format_mo_DP'] >= 8) c_Ash_QUAL_low = df_re_snp['QUAL'] < 30 c_Ash_snp_FS_bad = df_re_snp['info_FS'] >= 45 c_Ash_snp_FS_mid = (df_re_snp['info_FS'] >= 25) & (~c_Ash_snp_FS_bad) 
c_Ash_snp_QD_bad = df_re_snp['info_QD'] < 2.5
c_Ash_snp_QD_mid = (df_re_snp['info_QD'] < 4) & (~c_Ash_snp_QD_bad)

# filter SNPs df_re_snp.shape #df_temp_snp = df_re_snp[c_GT & c_AF & ~c_Ash_QUAL_low & ~c_Ash_snp_FS_bad & ~c_Ash_snp_FS_mid & ~c_Ash_snp_QD_bad & ~c_Ash_snp_QD_mid] df_temp_snp = df_re_snp[c_GT & c_AF & c_FILTER & c_DP]
df_temp_snp.shape
print list(df_temp_snp.columns)
df_temp_snp.ID.isnull().sum()
df_temp_snp[['family_id','ID', 'CHROM', 'POS', 'REF','ALT','info_ANN', 'score', 'format_p1_ref_AD', 'format_p1_alt_AD',
    'format_p1_DP','format_fa_DP', 'format_mo_DP']][df_temp_snp['score'] >= .75]

my_list = []
for sc in [0,.3, .4, .5, .6, .7, .75, .8]:
    print sc
    mysumm =  df_temp_snp.groupby(['fam']).apply(protsummNum, ('blood', 'saliva', 'prot',  sc))
    mysumm['score'] = sc
    my_list.append(mysumm)

len(my_list)


rerun_bl_v_sa_snp = pd.concat(my_list)
rerun_bl_v_sa_snp
rerun_bl_v_sa_snp.to_excel('/mnt/xfs1/home/asalomatov/nov19_rerun_bl_v_sa_snp.xls')
rerun_bl_v_sa_snp.to_excel('/mnt/xfs1/home/asalomatov/nov19_rerun_bl_v_sa_snp_noRS.xls')
rerun_bl_v_sa_snp.to_excel('/mnt/xfs1/home/asalomatov/nov19_rerun_bl_v_sa_snp_parents00.xls')
rerun_bl_v_sa_snp.to_excel('/mnt/xfs1/home/asalomatov/nov19_rerun_bl_v_sa_snp_parents00_noRS.xls')

### if desired remove snps with rs IDs
df_re_indel.shape
df_re_indel =  df_re_indel[df_re_indel.ID.isnull()]
df_re_indel.shape
#indel
c_GT = (df_re_indel['p1_gt_type'] >= 1) & (df_re_indel['fa_gt_type'] <= 0) & (df_re_indel['mo_gt_type'] <= 0)
c_AF = (df_re_indel['info_dbNSFP_1000Gp3_AF'] < .001) | (df_re_indel['info_dbNSFP_ExAC_AF'] < .001) | (df_re_indel['info_dbNSFP_ExAC_AF'].isnull())
c_FILTER = (df_re_indel['FILTER'] == 'PASS')
c_DP = (df_re_indel['format_p1_DP'] >= 8) & (df_re_indel['format_fa_DP'] >= 8) & (df_re_indel['format_mo_DP'] >= 8)
c_Ash_QUAL_low = df_re_indel['QUAL'] < 30
c_Ash_indel_QD_low = df_re_indel['info_QD'] < 1
sum(c_Ash_indel_QD_low)
c_Ash_indel_FS_bias = df_re_indel['info_FS'] >= 25
sum(c_Ash_indel_FS_bias)
c_Ash_indel_RP_bias = df_re_indel['info_ReadPosRankSum'] <= -3.
sum(c_Ash_indel_RP_bias)
# filter INDELs
df_re_indel.shape
df_temp_indel = df_re_indel[c_GT & c_AF & c_FILTER & c_DP & ~c_Ash_QUAL_low & ~c_Ash_indel_QD_low & ~c_Ash_indel_FS_bias & ~c_Ash_indel_RP_bias]
df_temp_indel['score'] = 0.0
c_p1_AB = (df_temp_indel['format_p1_alt_AD'].astype(float)/df_temp_indel['format_p1_DP'])  > .2
df_temp_indel.shape
df_temp_indel =df_temp_indel[c_p1_AB]
df_temp_indel.shape
print list(df_temp_indel.columns)
sf_indel_denovo = df_temp_indel[['CHROM', 'POS', 'REF','ALT','info_ANN', 'fam', 'prot', 'source']]
sf_indel_denovo['memb'] = 'p1'
sf_indel_denovo['source'] = 'sf'
sf_indel_denovo['descr'] = sf_indel_denovo['info_ANN'].apply(lambda x: x.split("|")[1])
sf_indel_denovo['gene'] = sf_indel_denovo['info_ANN'].apply(lambda x: x.split("|")[3])
sf_indel_denovo.drop('info_ANN', axis=1, inplace=True)
sf_indel_denovo
col_names = [i.lower() for i in sf_indel_denovo.columns]
col_names[0] = 'chr'
sf_indel_denovo.columns = col_names


df_temp_indel.groupby(['fam']).apply(protsummNum, ('bl', 'sa', 'prot', -1)).to_excel('/mnt/xfs1/home/asalomatov/nov19_rerun_bl_v_sa_indel.xls')
df_temp_indel.groupby(['fam']).apply(protsummNum, ('bl', 'sa', 'prot', -1)).to_excel('/mnt/xfs1/home/asalomatov/nov19_rerun_bl_v_sa_indel_noRS.xls')

df_bay_codified.columns = ['chr', 'pos', 'ref', 'alt', 'descr', 'gene', 'fam', 'memb', 'prot']



df_re_snp['status'][cAsh_snp_1 & cAsh_snp_2 & cAsh_snp_3].value_counts()
df_re_snp['dbsnp'][cAsh_snp_1 & cAsh_snp_2 & cAsh_snp_3].value_counts()
df_re_snp['dbsnp'].value_counts()


df_re_snp = pd.merge(df_re_snp, df_kn[['var', 'status', 'descr']], how='left', on='var')
df_re_snp.head()
df_re_snp['status'][df_re_snp.status.isnull()] = 'extra'
df_re_snp['dbsnp'] = 'yes'
df_re_snp['dbsnp'][df_re_snp.ID.isnull()] = 'no'
df_re_snp.ID.isnull().sum()



### NYG denovos
df_list = []
for f in [1,2,3,4]:
    for p in ['blood', 'saliva']:
        file_name = str(f)+'-'+p+'-denovo.txt'
        print file_name
        df = pd.read_csv("/mnt/scratch/asalomatov/BloodVsSaliva/"+file_name, header=None, sep="\t")
        df.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filt', 'info', 'fmt', 'gt_p1', 'gt_fa', 'gt_mo']
        df.head()
        df['gt_p1'] = df.gt_p1.apply(lambda df: df[:3])
        df['gt_fa'] = df.gt_fa.apply(lambda df: df[:3])
        df['gt_mo'] = df.gt_mo.apply(lambda df: df[:3])
        df['fam'] = str(f)
        df['prot'] = p
        df['memb'] = 'p1'
        df_list.append(df)

len(df_list)
for i in df_list:
    print i.shape

df_nyg = pd.concat(df_list)
df_nyg.head()
df_nyg['vartype'] = df_nyg.apply(varType, axis=1)
df_nyg.head()
df_nyg.dtypes
df_nyg.chr.value_counts()
df_nyg.gt_p1.value_counts()
df_nyg.gt_fa.value_counts()
df_nyg.gt_mo.value_counts()
df_nyg.isnull().sum()
df_nyg['pos_var'] = df_nyg['chr'].map(str) + '_' + df_nyg.pos.map(str)
df_nyg['fam_prot_pos_var'] = df_nyg['fam'] + '_' + df_nyg['prot'] + '_' + df_nyg['chr'].map(str) + '_' + df_nyg.pos.map(str)
df_nyg['source'] = 'nyg'
df_nyg['score'] = 0
df_nyg.shape
df_nyg = df_nyg[df_nyg.gt_p1.isin(['1/1','0/1']) &  df_nyg.gt_fa.isin(['0/0','0/0']) & df_nyg.gt_mo.isin(['0/0','0/0'])]
df_nyg.shape
df_nyg.head()

sum(501248 == df_nyg.pos)
'21854332'  in df_nyg.pos

df_nyg.groupby(['vartype', 'fam']).apply(protsummNum, ('blood', 'saliva', 'prot', -1))
df_nyg.groupby(['vartype', 'fam']).apply(protsummNum, ('blood', 'saliva', 'prot', -1)).to_excel('/mnt/xfs1/home/asalomatov/nov24_nyg_denovo_bl_v_sa.xls')

