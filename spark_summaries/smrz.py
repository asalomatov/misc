# this is for interactive work

def sumGene(x):
    Count = len(x)
    Count_LOF = x.lof.sum()
    Count_MIS = x.missense.sum()
    Count_SYN = x.syn.sum()
    Count_IMPACT_LOF = x.impact_lof.sum()
    Count_DMIS = (x.missense & x.dmg_miss).sum()
    N_indiv = round(len(x.SP_id.unique()), 0)
    indiv = '_'.join(x.SP_id.unique().tolist())
    SFARI_score = x.SFARIscore.iloc[0]
    lof_z_perc_rank = round(x.lof_z_perc_rank.iloc[0], 3)
    mis_z_perc_rank = round(x.mis_z_perc_rank.iloc[0], 3)
    pLI_perc_rank = round(x.pLI_perc_rank.iloc[0], 3)
    asd_score_perc_rank = round(x.asd_score_perc_rank.iloc[0], 3)
    LGDscore_perc_rank = round(x.LGDscore_perc_rank.iloc[0], 3)
    RVIS_perc_rank = round(x.RVIS_perc_rank.iloc[0], 3)
    return pandas.Series(
        [Count_LOF, Count_MIS, Count_SYN,
         Count_IMPACT_LOF, Count_DMIS, N_indiv,
         SFARI_score,
         lof_z_perc_rank, mis_z_perc_rank,
         pLI_perc_rank,
         asd_score_perc_rank,
         LGDscore_perc_rank,
         RVIS_perc_rank,
         indiv],
        index=['N_LOF', 'N_MIS', 'N_SYN',
               'N_IMPACT_LOF', 'N_DMIS',
               'N_indiv', 'SFARI',
               'lof_z', 'mis_z',
               'pLI',
               'asd(Olga)',
               'LGD',
               'RVIS',
               'indiv']
    )


def sumBatch(x):
    Count_trios = len(x.SP_id.unique())
    Count_LOF = sum(x.lof.astype(str) == 'True')
    Count_MIS = sum(x.missense.astype(str) == 'True')
    Count_SYN = sum(x.syn.astype(str) == 'True')
    Count_IMPACT_LOF = sum(x.impact_lof.astype(str) == 'True')
    Count_DMIS = sum(x.dmis.astype(str) == 'True')
    return pandas.Series(
        [Count_trios, Count_LOF, Count_MIS, Count_SYN,
         Count_IMPACT_LOF, Count_DMIS ],
        index=['N', 'N_LOF', 'N_MIS',
               'N_SYN', 'N_IMPACT_LOF', 'N_DMIS']
    )
# read denovo calls

import pandas
from natsort import natsorted

spark_dups = pandas.read_csv('/mnt/xfs1/scratch/asalomatov/data/SPARK/info/b1-10_dups.csv')
trio_dups = spark_dups[(spark_dups.batch.isin([4, 10]))]
ped.fam_id[ped.lab_id.isin(trio_dups.lab_id)].sort_values()
complete_dups_trio_id = ['A00062', 'A00073']
ped[ped.fam_id.isin(complete_dups_trio_id)]
complete_dups_sp_id = ['SP0001293', 'SP0001681']

ped_last = pandas.read_csv('/mnt/xfs1/scratch/asalomatov/data/SPARK/info/b10_pedigree.csv')
ped_last['batch'] = 10
ped_old = pandas.read_table('/mnt/xfs1/scratch/asalomatov/data/SPARK/info/b1-9_pedigree.tsv')
ped_old.shape[0]/3
ped_old.columns
ped = pandas.concat([ped_old[ped_last.columns], ped_last ])
ped.batch = 'b' + ped.batch.astype(str)
ped.ix[ped.batch == 'b1', 'batch'] = 'b1-2'
ped.ix[ped.batch == 'b2', 'batch'] = 'b1-2'
ped = ped[~((ped.batch == 'b4') & ped.fam_id.isin(complete_dups_trio_id))]

n_trios_by_b = ped.groupby('batch').fam_id.apply(lambda i: len(set(i)))
print n_trios_by_b
print n_trios_by_b.sum()


dnv = pandas.read_csv('/mnt/ceph/users/asalomatov/spark/denovo/def1/SF_denovo_b1-2_b3_b4_b5_b6_b7_b8_b9_b10_all_research.csv')
# dnv = pandas.read_csv('/mnt/ceph/users/asalomatov/spark/denovo/results/b1-9/SF_denovo_b1-2_b3_b4_b5_b6_b7_b8_b9_all_research.csv')
dnv.head()
dnv.shape
dnv.dtypes
dnv.lof.value_counts().sum()
dnv.dmg_miss.value_counts()
dnv.missense.value_counts()
dnv.syn.value_counts()
dnv.batch.value_counts()
dnv['dmis'] = dnv.missense & dnv.dmg_miss
#now remove duplicated trios from batch 4
c_b = dnv.batch == 'b4'
c_tr = dnv.SP_id.isin(complete_dups_sp_id)
sum(c_tr)
sum(c_b & c_tr)
dnv[~(c_b & c_tr)].shape
dnv = dnv[~(c_b & c_tr)]

# summary by batch
def summarizeByBatch(df, trio_counts_by_batch=n_trios_by_b):
    dnv_by_batch = df.groupby(by='batch').apply(sumBatch)
    dnv_by_batch = dnv_by_batch.reindex(index=natsorted(dnv_by_batch.index))
    dnv_by_batch = dnv_by_batch.assign(N_trios=n_trios_by_b)
    dnv_by_batch = dnv_by_batch[['N_trios', 'N_LOF', 'N_MIS', 'N_SYN', 'N_IMPACT_LOF', 'N_DMIS']]
    return dnv_by_batch

dnv_by_batch = summarizeByBatch(dnv)

dnv_by_batch.to_csv('~/spark_sum_b1-10/vars_by_batch.csv')

dnv_by_batch.N_LOF.sum()
dnv_by_batch.sum()

# by batch in spark genes
dnv.shape
dnv.spark_genes.value_counts()
dnv.spark_genes.dtypes
c_spark = dnv.spark_genes
dnv_by_batch_clinical = summarizeByBatch(dnv[c_spark])
dnv_by_batch_clinical
dnv_by_batch_clinical.sum()
dnv_by_batch_clinical.to_csv('~/spark_sum_b1-10/vars_by_batch_clinical.csv')

# summary by gene
dnv_by_gene = dnv.groupby('ANN.GENE').apply(lambda row: sumGene(row))
dnv_by_gene.head()
dnv_by_gene['N_LOFMIS'] = dnv_by_gene.N_LOF + dnv_by_gene.N_MIS
dnv_by_gene.sort_values('N_LOFMIS', inplace=True, ascending=False)
dnv_by_gene.head()

dnv_by_gene.shape

dnv_by_gene.to_csv('~/spark_sum_b1-10/lofmis_by_gene.csv')

# summary by gene in spark genes
dnv_by_gene_clinical = dnv[c_spark].groupby('ANN.GENE').apply(lambda row: sumGene(row))
dnv_by_gene_clinical['N_LOFMIS'] = dnv_by_gene_clinical.N_LOF + dnv_by_gene_clinical.N_MIS
dnv_by_gene_clinical.sort_values('N_LOFMIS', inplace=True, ascending=False)
dnv_by_gene_clinical.head()

dnv_by_gene_clinical.shape

dnv_by_gene.to_csv('~/spark_sum_b1-10/lofmis_by_gene_clinical.csv')

