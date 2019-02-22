import sys
import gzip

with open(sys.argv[3], 'r') as blkf:
	blkd={}
	for line in blkf:
		blkd[line[:-1]]=None

def opener(filename):
    f = open(filename,'rb')
    if (f.read(2) == '\x1f\x8b'):
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f

inp=opener(sys.argv[1])



log=open('.'.join(sys.argv[1].split('.')[:-1])+'.filter_reads_out.log', 'w')

if sys.argv[2]=='-': 
	outp=sys.stdout
	outpname='stdout'
else:
	outpname=sys.argv[2] 
	outp=open(outpname, 'w') 
log.write('Filtering reads (matching those listed in '+sys.argv[3]+') out from '+sys.argv[1]+'\n')
cnttot=cntrem=cntout=0
for line in inp:
	cnttot=cnttot+1
	if line.split()[0] not in blkd: outp.write(line); cntout=cntout+1
cntrem=cnttot-cntout
log.write('# reads in '+sys.argv[1]+': '+str(cnttot)+'\n# reads in '+sys.argv[3]+': '+str(len(blkd))+'\n# reads removed: '+str(cntrem)+'\n# reads left after filtering and written to '+outpname+': '+str(cntout)+'\n')
inp.close()
log.close()
outp.close()

