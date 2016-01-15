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

temp_list_snp_nyg = []
temp_list_indel_nyg = []
for fam in ['1', '2', '4']:
    for prot in ['blood', 'saliva']:
        fam_id = fam + '-' + prot
        print fam_id
        myfile = '/mnt/scratch/asalomatov/BloodVsSaliva/'+fam_id+'.snpeff.vcftools-target.vcf.gz'
        print  myfile
        myvars = variants.Variants(myfile, fam_id)
        myvars.readFromVcf()
        df_temp =  myvars.variants
        #df_temp.columns = [i.replace(fam+'.', '') for i in df_temp.columns]
        #print df_temp.columns
        df_temp['pos_var'] = df_temp['CHROM'].map(str) + '_' + df_temp.POS.map(str)
        df_temp['fam'] = df_temp['family_id'].apply(lambda x: x.split('-')[0])
        df_temp['prot'] = df_temp['family_id'].apply(lambda x: x.split('-')[1])
        df_temp['source'] = 'rerun'
        #df_temp = df_temp[df_temp.ID.isnull()]
        temp_list_indel_nyg.append(df_temp[~df_temp.vartype.isin(['snp'])])
        temp_list_snp_nyg.append(df_temp[~df_temp.score.isnull() & df_temp.vartype.isin(['snp'])])

