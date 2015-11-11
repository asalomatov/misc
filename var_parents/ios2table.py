#!/usr/bin/python
'''
look up sequence info for SSCexome.txt  coordinates using pysam, extract some columns
'''
import pysam, sys

#if len(sys.argv) != 4:
#    print 'usage:'
#    print sys.argv[0], '/pathto.ref.fa /path/to/input.bed /path/to/output.vcf'
#    sys.exit(1)
#
#genref, iosfile = sys.argv[1:]
#genref = sys.argv[1]
iosfile = '../../data/SSCdeNovoCalls/SSCexome.csv'
genref = '/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'

def transformIosLine(x, ref = genref):
    '''
    x is a list with a tab split line from SSCexome.txt file. Transform it to a minimal vcf line
    retaining all of the information in an infor field.
    '''
    fam_id = x[0]
    study = x[1]
    chrom, ios_coord = x[2].split(':')
    vcf_coord = ios_coord
    variant = x[3]
    #misc_info = x[4:]
    #misc_info.append(study)
    #misc_info = '__'.join(map('_'.join, misc_info))
    #misc_info = '__'.join(['_'.join(x) for x in  misc_info])
    var_desc = variant[0:3]
    ios_ref = ''
    ios_alt = ''
    if var_desc == 'sub':
        y = variant.split('->') 
        ios_ref = y[0][-1]
        ios_alt = y[1][0]
    elif var_desc == 'ins':
        vcf_coord = str(int(ios_coord) - 1)

    vcf_id = '.'
    qual = '.'
    filt = '.'
    allels = extracAllels(chrom, vcf_coord, variant, genref)
    geno = x[4] 
    gn_ref, gn_alt = geno.split('/')
    gn_ref = gn_ref.split()
    gn_alt = gn_alt.split()
#    print gn_ref, gn_alt
    gn_parents = int(gn_alt[0]) + int(gn_alt[1])
    if gn_parents < 3:
        return ''
    gn_p1 = int(gn_alt[2])
    gn_s1 = 'na'
    if len(gn_ref) == 3: #trio
        if (gn_parents == 3 and gn_p1 == 0) or (gn_parents == 4 and gn_p1 <= 1):
            return '\t'.join([fam_id, chrom, vcf_coord, allels ])
        else: 
            return ''
    else:
        gn_s1 = int(gn_alt[3])
        gn_kids = gn_p1 + gn_s1
        if (gn_parents == 3 and gn_kids <= 1) or (gn_parents == 4 and gn_kids <= 2):
            return '\t'.join([fam_id, chrom, vcf_coord, allels ])
        else: 
            return ''
#    else: # quad
#        return '\t'.join([fam_id, chrom, vcf_coord, allels])
#    return '\t'.join([chrom, vcf_coord, vcf_id, allels, qual, filt, 'FAM='+fam_id+';'])
#    info = '.'
#    return '\t'.join([chrom, vcf_coord, vcf_id, allels, qual, filt, info])
    
def extracAllels(chrom, vcf_coord, var_descr, genref):
    '''compute alternative allel from description in bed file.
    must be one of sub(C->T), ins(CCCT), del(5)
    '''
    ref_allel = pysam.faidx(genref, chrom+':'+vcf_coord+'-'+vcf_coord)[1].strip()
    if 'sub' in var_descr:
        yy = var_descr.split('->')
        ref = yy[0][-1]
        if ref == ref_allel:
            return '\t'.join([ref, yy[1][0]])
        else:
            print 'ref allels do not match, exiting...'
            print chrom, vcf_coord, var_descr
            sys.exit(1)
    elif 'ins' in var_descr:
        yy = var_descr.split('(')[1]
        return '\t'.join([ref_allel, ref_allel + yy[:-1]])
    elif 'del' in var_descr:
        yy = var_descr.split('(')[1]
        del_len = int(yy[:-1])
        vcf_coord_end = str(int(vcf_coord) + int(del_len))
        ref = pysam.faidx(genref, chrom+':'+vcf_coord+'-'+vcf_coord_end)[1].strip()
        if ref[0] == ref_allel:
            return '\t'.join([ref, ref_allel])
        else:
            print 'ref allels do not match in del, exiting...'
            print chrom, vcf_coord, var_descr
            sys.exit(1)
    else:
        print 'format not found, exiting...'
        sys.exit(1)

#for l in sys.stdin:
#    l_list = l.split('\t')
#    misc_info = ['_'.join(x.split()) for x in l_list]
#    misc_info = '__'.join(misc_info)
#    misc_info = misc_info.replace(';','___',100)
#    vcfline = transformIosLine(l_list)
#    print vcfline+'MISC_INFO='+misc_info

with open(iosfile, 'r') as ios:
    for i, l in enumerate(ios):
        if i < 2:
            continue
        l_list = l.split('\t')
        misc_info = ['_'.join(x.split()) for x in l_list]
        misc_info = '__'.join(misc_info)
        misc_info = misc_info.replace(';','___',100)
        outputline = transformIosLine(l_list)
        if outputline:
            print '\t'.join([outputline, misc_info])
#        if i > 10 : sys.exit(1)
