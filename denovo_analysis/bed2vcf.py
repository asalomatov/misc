'''
look up sequence info for bed coordinates using pysam, write as headerless vcf.
'''
import pysam, sys

if len(sys.argv) != 4:
    print 'usage:'
    print sys.argv[0], '/pathto.ref.fa /path/to/input.bed /path/to/output.vcf'
    sys.exit(1)

genref, bedfile, outvcf = sys.argv[1:]

def transformBedLine(x, ref = genref):
    chrom = x[0]
    vcf_id = '.'
    vcf_pos = str(int(x[1]) + 1)
    ref_allel = pysam.faidx(genref, x[0]+':'+vcf_pos+'-'+vcf_pos)[1].strip()
    allels = extracAllels(chrom, vcf_pos, ref_allel, x[3], genref)
    qual = '.'
    filt = '.'
    info = '.'
    return '\t'.join([chrom, vcf_pos, vcf_id, allels, qual, filt, info])
    
def extracAllels(chrom, vcf_pos, ref_allel, var_descr, genref):
    '''compute alternative allel from description in bed file.
    must be one of sub(C->T), ins(CCCT), del(5)
    '''
    if 'sub' in var_descr:
        yy = var_descr.split('->')
        ref = yy[0][-1]
        if ref == ref_allel:
            return '\t'.join([ref, yy[1][0]])
        else:
            print 'ref allels do not match, exiting...'
            print chrom, vcf_pos
            sys.exit(1)
    elif 'ins' in var_descr:
        yy = var_descr.split('(')[1]
        return '\t'.join([ref_allel, ref_allel + yy[:-1]])
    elif 'del' in var_descr:
        yy = var_descr.split('(')[1]
        del_len = int(yy[:-1])
        vcf_pos_end = str(int(vcf_pos) + int(del_len))
        ref = pysam.faidx(genref, chrom+':'+vcf_pos+'-'+vcf_pos_end)[1].strip()
        if ref[0] == ref_allel:
            return '\t'.join([ref, ref_allel])
        else:
            print 'ref allels do not match in del, exiting...'
            print crom, vcf_pos
            sys.exit(1)
    else:
        print 'format not found, exiting...'
        sys.exit(1)


with open(outvcf, 'a') as fout:
    with open(bedfile, 'r') as bed:
        for i, l in enumerate(bed):
            l_list = l.split()
            vcfline = '\t'.join([transformBedLine(l_list), l_list[-1]])
#            print vcfline
#            if i > 10 : sys.exit(1)
            fout.write(vcfline+'\n')


