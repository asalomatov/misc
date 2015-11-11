import sys

num_mism = 0
for i, l in enumerate(sys.stdin):
    l_list = l.rstrip('\n').split('\t')
    if i > 0:
        hg38_chr = l_list[0]
        hg38_pos = l_list[1]
        hg19_chr = l_list[7]
        hg19_pos = l_list[8]
        if hg38_chr + hg38_pos != hg19_chr + hg19_pos:
            num_mism += 1
        l_list[0] = hg19_chr
        l_list[1] = hg19_pos
        l_list[7] = hg38_chr
        l_list[8] = hg38_pos
    print '\t'.join(l_list)
    if i % 1000000 == 0:
        print >> sys.stderr, 'line number ', i
print >> sys.stderr, 'number of mismatches is ', num_mism
