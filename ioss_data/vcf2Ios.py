#!/usr/bin/python
'''
Convert a vcf file to txt for comparison with Iosifov's data
'''
import sys, os

if len(sys.argv) != 2:
    print 'Usage:'
    print sys.argv[0], '/path/to/input.vcf outputdir'


inf, outdir = sys.argv[1:]
print 'input file: ', inf
outf = os.path.join(outdir, os.path.basename(inf)+'.txt')
print 'output file: ', outf

with open(outf, 'w') as fout:
    with open(inf, 'r') as fin:
        for l in fin:
            if l[0] == '#':
                continue
            ls = l.split()
            print ls
            sys.exit(1)
            if ls[0] == " ":
                fout.write('\t'.join(ls)+'\n')


