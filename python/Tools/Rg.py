import re
import sys
import argparse

parser = argparse.ArgumentParser(prog='Rg-fromProp',description='Compute average gyration radius from a trjProp output')
parser.add_argument('--input','-i', action = 'store', type = str, help = 'Input .gyr file', required=True)
parser.add_argument('--begin','-b', action = 'store', type = float, help = 'Where to begin to make average ', required=False,default=0.0)

args = parser.parse_args()
begin=args.begin;
try:
    fp=open(args.input,'r')
except FileNotFoundError:
    print(args.input+": File not found!")
    sys.exit(2)

findTot=re.compile(r'^Tot')
findPol=re.compile(r'^Pol')
findNoH=re.compile(r'^NoH')
findClust=re.compile(r'^.* (\d\d*)  *time')
findRg=re.compile(r'^.* Rg = *(.*)  *a =')
findTime=re.compile(r'^.* time = *(.*)  *Rg =')
lines=fp.readlines()
Rg_t={}
Rg_p={}
Time=[]
for line in lines:
    if findTot.match(line):
        Lab=findClust.match(line).group(1)
        Rg=float(findRg.match(line).group(1))
        Time.append(float(findTime.match(line).group(1)))
        if len(Rg_t) == 0:
            Rg_t[Lab]=[]
        Rg_t[Lab].append(Rg)
    elif findPol.match(line) or findNoH.match(line):
        Lab=findClust.match(line).group(1)
        Rg=float(findRg.match(line).group(1))
        if len(Rg_p) == 0:
            Rg_p[Lab]=[]
        Rg_p[Lab].append(Rg)
#   Average over Tot
print(len(Time))
Tmax=len(Time)-len([i for i in Time if i > begin])
print(Tmax)
RgT={}
RgP={}
for Rg in Rg_t:
    RgT[Rg]=sum(Rg_t[Rg][Tmax:])
    RgT[Rg]/=len(Rg_t[Rg][Tmax:])



for Rg in Rg_p:
    RgP[Rg]=sum(Rg_p[Rg][Tmax:])
    RgP[Rg]/=len(Rg_p[Rg][Tmax:])

for Rg in RgT:
    print(" Rg Tot = %12.3f " % RgT[Rg])
for Rg in RgP:
    print(" Rg NoH = %12.3f " % RgP[Rg])
