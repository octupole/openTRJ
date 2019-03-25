import re
import sys
import argparse

parser = argparse.ArgumentParser(prog='Rg-fromProp',description='Compute average gyration radius from a trjProp output')
parser.add_argument('--input','-i', action = 'store', type = str, help = 'Input .gyr file', required=True)

args = parser.parse_args()
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
        if len(Rg_t) == 0:
            Rg_t[Lab]=[]
        Rg_t[Lab].append(Rg)
        Time.append(float(findTime.match(line).group(1)))
    elif findPol.match(line) or findNoH.match(line):
        Lab=findClust.match(line).group(1)
        Rg=float(findRg.match(line).group(1))
        if len(Rg_p) == 0:
            Rg_p[Lab]=[]
        Rg_p[Lab].append(Rg)
#   Average over Tot

#   Average over Tot
RgT={}
for Rg in Rg_t:
   for o in range(len(Rg_t[Rg])):
       print(" %10.1f %10.3f  %10.3f "%(Time[o],Rg_t[Rg][o],Rg_p[Rg][o]))
