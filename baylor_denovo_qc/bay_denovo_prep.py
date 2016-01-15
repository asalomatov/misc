import pandas as pd
import numpy as np
import os

bay_info = pd.read_csv("/mnt/scratch/asalomatov/baylor_qc/Simons_pilot_sample_info", header=False, sep="\t")
bay_info.head()
def myf(row):
    if row['Grant ID'][0]=='A':
        return 3
    elif row['Grant ID'][0]=='P':
        return 1
    elif row['Grant ID'][0]=='M':
        return 2
    else:
        return 4

def myf1(row):
    if row['Relationship']=='Proband':
        return 'p1'
    elif row['Relationship']=='Father':
        return 'fa'
    elif row['Relationship']=='Mother':
        return 'mo'
    elif row['Relationship']=='Son1':
        return 'p1'
    elif row['Relationship']=='Son2':
        return 's1'
    elif row['Relationship']=='MatGrandMom':
        return 'gm'
    elif row['Relationship']=='MatGrandDad':
        return 'gf'
    else:
        return None

def myf2(row):
    if row['Sample Type']=='Blood':
        return 'bl'
    elif row['Sample Type']=='Saliva':
        return 'sa'
    else:
        return None

def myf3(row):
    if row['memb']=='p1' or row['memb']=='s1':
        return str(row['fam_id'])+'-fa'
    else:
        return '0'

def myf4(row):
    if row['memb']=='p1' or row['memb']=='s1':
        return str(row['fam_id'])+'-mo'
    else:
        return '0'


bay_info.Relationship.value_counts()
bay_info.apply(myf, axis=1)
bay_info['fam'] = bay_info.apply(myf, axis=1)
bay_info.fam.value_counts()
bay_info['memb'] = bay_info.apply(myf1, axis=1)
bay_info.memb.value_counts()
bay_info['prot'] = bay_info.apply(myf2, axis=1)
bay_info.prot.value_counts()
bay_info.head()
bay_info['fam_id'] =  bay_info['fam'].astype(str) + '-' + bay_info['prot']
bay_info['smpl_id'] =  bay_info['fam_id'] + '-' + bay_info['memb']
bay_info['father'] = bay_info.apply(myf3, axis=1)
bay_info['mother'] = bay_info.apply(myf4, axis=1)
bay_info['sex'] = 1
bay_info['sex'][bay_info['memb']=='mo'] = 2
bay_info['sex'][bay_info['memb']=='gm'] = 2
bay_info['sex'][(bay_info['memb']=='p1') & (bay_info['fam']==1)] = 2
bay_info.head()
bay_info.sex.value_counts()
bay_info['pheno'] = 1
bay_info['pheno'][bay_info['memb']=='p1'] = 2
bay_info.pheno.value_counts()
bay_info['bam'] = bay_info.apply(lambda row: '/'.join(['/mnt/scratch/asalomatov/baylor_qc', row['BAM file name'].split('/')[1]]), axis=1)
bay_info.bam.value_counts()
bay_info['bay_smpl_id'] =  bay_info['Grant ID'] + '_' + bay_info['Relationship'] + '_' + bay_info['Sample Type']
bay_info.bay_smpl_id.value_counts()
bay_info.head()
bay_info.isnull().sum()
# sf id to bay_id dict
sf_bay = {}
for i, row in bay_info.iterrows():
    sf_bay[row['smpl_id']] = row['bay_smpl_id']

# write ped
bay_info[['fam_id','smpl_id','father','mother','sex','pheno']].to_csv('/mnt/scratch/asalomatov/baylor_qc/dnm_filter/baylor.ped', sep="\t",header=False, index=False)
# write dnm bam file map
bay_info[['smpl_id','bam']].to_csv('/mnt/scratch/asalomatov/baylor_qc/dnm_filter/baylor-bams.csv', sep="\t", header=False, index=False)

# split baylor vcf by family
import commands
work_dir = '/mnt/scratch/asalomatov/baylor_qc/'
for f in [1,2,3,4]:
    for prot in ['bl', 'sa']:
        if f==3 and prot=='sa':
            continue
        ofile = os.path.join(*[work_dir,'dnm_filter',str(f)+'-'+prot+'.txt'])
        mysamples = ','.join([sf_bay[str(f)+'-'+prot+'-p1'],sf_bay[str(f)+'-'+prot+'-fa'],sf_bay[str(f)+'-'+prot+'-mo']])
        cmd = 'zcat '+ os.path.join(work_dir,'simons_pvcf_snps.annotated.onTarget.vcf.gz')+ ' | vcf-subset -r -e -c '+mysamples+' | grep -v ^# > '+ofile
        #cmd = 'zcat '+ os.path.join(work_dir,'simons_pvcf_snps.annotated.onTarget.vcf.gz')+ ' | vcf-subset -r -e -c '+mysamples+' > '+ofile
        print commands.getstatusoutput(cmd)




