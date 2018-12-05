import sys
import re

class xvgFile:
    def __init__(self, args,cutoff=100):
        self.cut=cutoff
        myxvg=''.join(args)
        self.xvgs=re.sub(r'@.*\n','',myxvg).split('&')
        self.xvgs.sort(key=lambda a: len(a),reverse=True)

    def __repr__(self):
        outf=''
        M=0
        for id0 in self.xvgs:

           if len(id0) != 1 and len(id0.split('\n')) > self.cut:    
                outf+='@    s%-d legend "Id No. %-s"\n' % (M,M)
                outf+='@    s%-d symbol 1 \n' % (M)
                outf+='@    s%-d symbol size 0.35 \n' % (M)
                outf+='@    s%-d symbol fill pattern 1\n' % (M)
                outf+='@    s%-d line type 0 \n' % (M)
                M+=1
        for id0 in self.xvgs:
           if len(id0) != 1 and len(id0.split('\n')) > self.cut:
                outf+=id0+'&'
        return outf
    



files=['_Ag.out','_Cl.out','_E0.out','_E1.out','_Rg.out']
files_n=[]
if len(sys.argv) == 1:
    print('No argument given abort')
    sys.exit()

program_name=sys.argv[0]
label=sys.argv[1]
cutoff=100
if len(sys.argv) > 2:
    cutoff=int(sys.argv[2])

for n in range(len(files)):
    files_n.append(label+'N'+files[n])
    files[n]=label+files[n]

infiles=[]
ofiles=[]
for file in files_n:
    ofiles.append(open(file,'w'))

for file in files:
    try:
        fn=open(file,'r')
        infiles.append(fn)
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.exit()

for n in range(len(infiles)):
    infile=infiles[n].readlines()
    outstr=xvgFile(infile,cutoff)
    ofiles[n].write(repr(outstr))
    sys.exit()
