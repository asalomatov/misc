#!/usr/bin/python

import sys

#infile, outfile = sys.argv[1:]

for l in sys.stdin:
    x = l.rstrip('\n').split('\t')
    fam_id = x[0]
    chrom, pos, ref, alt = x[3].split(':')
    child = x[4]
    misc_info = '__'.join(x)
    vcf_id = '.'
    qual = '.'
    filt = '.'
    vcf_line = '\t'.join([chrom, pos, vcf_id, ref, alt, qual, filt, 'FAM='+fam_id+';CHILD='+child+';MISC_INFO='+misc_info])
    print vcf_line
