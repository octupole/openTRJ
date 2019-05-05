
import re
import sys
import argparse
import math
import numpy as np

Hydrophobic={'ARG':'C','LYS':'C','ASP':'C','GLU':'C',
'GLN':'P','ASN':'P','HIS':'P','HSD':'P','SER':'P','THR':'P','TYR':'P','CYS':'P',
'ALA':'H','ILE':'H','LEU':'H','MET':'H','PHE':'H','VAL':'H','PRO':'H','GLY':'H',
'SOL':'SOL','NA':'NA'
}

parser = argparse.ArgumentParser(prog='RadialAvg',description='Average radial density')
parser.add_argument('--input','-i', action = 'store', type = str, help = 'Input .rho file', required=True)
parser.add_argument('--output','-o', action = 'store', type = str, help = 'output .xvg file', 
required=False,default='output.xvg')
parser.add_argument('--smooth','-s', action = 'store', type = int, help = 'Smooth data', 
required=False,default=1)

args = parser.parse_args()

N=args.smooth
if N == 0:
    N=1

filein=args.input
fileout=args.output
try:
    fp=open(filein,'r')
except FileNotFoundError:
    print(args.input+": File not found!")
    sys.exit(2)
fout=open(fileout,'w')

data={}
dataP={}
lines=fp.readlines()
nlegend=0
legend={}
Rad=0
mymod='valid'

Types={}
fout.write('@    s0 legend  "Total"\n')
fout.write('@    s1 legend  "Prot"\n')
M=2
for line in lines:
    if re.search('^@ *',line):
        amino=line.split()[-1].replace('"','')
        legend[nlegend]=amino
        Types[Hydrophobic[amino]]=1
        nlegend+=1

dataa={}
for labels in sorted(Types.keys()):
    print(labels)
    dataa[labels]={}
    fout.write('@    s%-2d legend  "%s"\n' % (M,labels))
    M+=1
print('Next should be identical to above')
for line in lines:
    if re.search('^&.*',line):
        Rad+=1
        continue
    if not re.search('^@  *',line):
        my=line.split()
        my1=float(my[0])
        my2=float(my[1])
        A=Hydrophobic[legend[Rad]]
        if my1 in dataa[A]:
            dataa[A][my1]+=my2
        else:
            dataa[A][my1]=my2
for A in sorted(dataa):
    for my in dataa[A]:
        if my in data:
            data[my]+=dataa[A][my]
        else:
            data[my]=dataa[A][my]
        if A == 'C' or A == 'P' or A == 'H':
            if my in dataP:
                dataP[my]+=dataa[A][my]
            else:
                dataP[my]=dataa[A][my]
        
x=[]
for d in sorted(data.keys()):
    x.append(data[d])
x=np.convolve(x,np.ones((N,))/N, mode=mymod)
x=np.append(x,np.zeros((N-1,)))
m=0
for d in sorted(data.keys()):
    fout.write(" %9.3f  %12.6f\n" %(d,x[m]))
    m+=1

fout.write('&\n')
x=[]
for d in sorted(dataP.keys()):
    x.append(dataP[d])
x=np.convolve(x,np.ones((N,))/N, mode=mymod)
x=np.append(x,np.zeros((N-1,)))
m=0
for d in sorted(data.keys()):
    fout.write(" %9.3f  %12.6f\n" %(d,x[m]))
    m+=1

fout.write('&\n')

for type in sorted(dataa.keys()):
    print(type)
    x=[]
    for d in sorted(dataa[type].keys()):
        x.append(dataa[type][d])
    x=np.convolve(x,np.ones((N,))/N, mode=mymod)
    x=np.append(x,np.zeros((N-1,)))
    m=0
    for d in sorted(dataa[type].keys()):
        fout.write(" %9.3f  %12.6f\n" %(d,x[m]))
        m+=1
    fout.write('&\n')
