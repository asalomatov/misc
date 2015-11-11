'''

'''
import sys, subprocess, commands
sys.path.insert(0, '/nethome/asalomatov/projects/pipeline/ppln/')
import logProc

print '\nsys.args   :', sys.argv[1:]
refGenome, tmpdir, gatk, dbsnp, outdir, outfile, inbed = sys.argv[1:8]

options = '''  \
--annotation BaseQualityRankSumTest  \
--annotation FisherStrand  \
--annotation GCContent  \
--annotation HaplotypeScore  \
--annotation HomopolymerRun  \
--annotation MappingQualityRankSumTest  \
--annotation MappingQualityZero  \
--annotation QualByDepth  \
--annotation ReadPosRankSumTest  \
--annotation RMSMappingQuality  \
--annotation DepthPerAlleleBySample  \
--annotation Coverage  \
--annotation ClippingRankSumTest  \
--annotation DepthPerSampleHC  \
-nct 20  
'''

I = ' -I '
inbams = ''
for f in sys.argv[8:]:
    inbams += I + f
print inbams
    
cmd = 'java -Xms750m -Xmx250g -XX:+UseSerialGC -Djava.io.tmpdir=%(tmpdir)s -jar %(gatk)s -T HaplotypeCaller %(inbams)s -o %(outfile)s -R %(refGenome)s --dbsnp %(dbsnp)s -L %(inbed)s %(options)s'
cmd = cmd % locals()
print cmd
logProc.logProc(outfile, outdir, cmd, 'started')
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
if p.returncode == 0:
    logProc.logProc(outfile, outdir, cmd, 'finished')
else:
    logProc.logProc(outfile, outdir, cmd, 'failed', stderr)

