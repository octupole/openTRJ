#!/Users/marchi/anaconda3/bin/python3
import re
import sys
import argparse
import sys

mass={'H': 1.008,'C': 12.0107,'O':15.9994,'N':14.0067,'NA':22.9897,'P':30.9738,'S':32.065,'CL':35.453,'K':39.0983,'CA':40.078}

parser = argparse.ArgumentParser(prog='molMass',description='Compute molecular mass from a pdb file')
parser.add_argument('--input0','-i0', action = 'store', type = str, help = 'Input xvg file', required=True)
parser.add_argument('--input1','-i1', action = 'store', type = str, help = 'Input xvg file', required=True)
parser.add_argument('--output','-o', action = 'store', type = str, help = 'Output avg xvg file', required=True)

args = parser.parse_args()
try:
    fp0=open(args.input0,'r')
except FileNotFoundError:
    print(args.input0+": File not found!")
    sys.exit(2)
try:
    fp1=open(args.input1,'r')
except FileNotFoundError:
    print(args.input1+": File not found!")
    sys.exit(2)

fout=open(args.output,'w') 

findPar=re.compile(r'^@  .*$')
findAmpers=re.compile(r'^&')

lines0=fp0.readlines()
lines1=fp1.readlines()

Rgs={}
for line in lines0:
    if not findPar.match(line) and not findAmpers.match(line):
        xy=line.split()
        x=int(float(xy[0]))
        y=float(xy[1])
        if x in Rgs.keys():
            Rgs[x].append(y)
        else:
            Rgs[x]=[]

Ags={}
for line in lines1:
    if not findPar.match(line) and not findAmpers.match(line):
        xy=line.split()
        x=int(float(xy[0]))
        y=float(xy[1])
        if x in Ags.keys():
            Ags[x].append(y)
        else:
            Ags[x]=[]

Rgs=sorted([[key,Rgs[key]] for key in Rgs.keys()])
Ags=sorted([[key,Ags[key]] for key in Ags.keys()])

for key in range(len(Rgs)):
    zRg=0.0
    Rg0=0.0
    agg=0
    zNorm=0.0
    Norm=0.0
    for o in range(len(Rgs[key][1])):
        Rg=Rgs[key][1][o]
        Agg=Ags[key][1][o]
        v=Rg*Rg*Rg
        zRg+=Rg*v
        zNorm+=v
        Rg0+=Rg
        Norm+=1.0
        agg+=Agg*v

    fout.write(" %d  %10.2f %12.5f %10.2f %12.5f \n" % (key,Norm, Rg0/Norm, agg/zNorm, zRg/zNorm,))
    