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

temp_list_snp = []
temp_list_indel = []
for fam in ['1', '2', '3', '4']:
    for prot in ['blood', 'saliva']:
        if prot == 'saliva' and fam == '3': continue
        fam_id = ''
        if 'bl' in prot:
            fam_id = fam+'-'+'bl'
        if 'sa' in prot:
            fam_id = fam+'-'+'sa'
        print fam_id
        #myfile = '/mnt/scratch/asalomatov/baylor_qc/rerun/'+prot+'/'+fam+'-HC-vars-flr-ann.vcf.gz'
        myfile = '/mnt/scratch/asalomatov/baylor_qc/rerun/'+prot+'/targ/'+fam+'-HC-ann-targ.vcf.gz'
        my_dnmfile = '/mnt/scratch/asalomatov/baylor_qc/rerun/dnm_filter/output/'+fam_id+'-dnm.csv'
        print  myfile
        print  my_dnmfile
        myvars = variants.Variants(myfile, fam_id)
        myvars.readFromVcf()
        df_temp =  myvars.variants
        df_temp.columns = [i.replace(fam+'.', '') for i in df_temp.columns]
        print df_temp.columns
        df_temp['pos_var'] = df_temp['CHROM'].map(str) + '_' + df_temp.POS.map(str)
        df_temp['fam'] = fam
        df_temp['prot'] = prot
        df_temp['source'] = 'SF'
        #df_temp = df_temp[df_temp.ID.isnull()]
        temp_list_indel.append(df_temp[~df_temp.vartype.isin(['snp'])])
        dnm_scores = pd.read_csv(my_dnmfile, header=None) 
        dnm_scores.columns = ['desc', 'chr', 'pos', 'score']
        dnm_scores['pos_var'] = dnm_scores['chr'].map(str) + '_' + dnm_scores.pos.map(str)
        df_temp = pd.merge(df_temp, dnm_scores[['pos_var', 'score']], how='left', on='pos_var')
        df_temp = df_temp[~df_temp.score.isnull() & df_temp.vartype.isin(['snp'])]
        temp_list_snp.append(df_temp)

