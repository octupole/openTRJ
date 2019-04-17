import re
import sys
import argparse

parser = argparse.ArgumentParser(prog='Rg-fromProp',description='Compute average gyration radius from a trjProp output')
parser.add_argument('--input','-i', action = 'store', type = str, help = 'Input .gyr file', required=True)
parser.add_argument('--begin','-b', action = 'store', type = float, help = 'Where to begin ', required=False,default=0.0)

args = parser.parse_args()
begin=args.begin;
try:
    fp=open(args.input,'r')
except FileNotFoundError:
    print(args.input+": File not found!")
    sys.exit(2)

findTot=re.compile(r'^Tot')
findClust=re.compile(r'^.* (\d\d*)  *time')
findTime=re.compile(r'^.* time = *(.*)  *Rg =')
lines=fp.readlines()
Time=[]
NCl=[]
newTime=-1
oldTime=-1
for line in lines:
    if findTot.match(line):
        newTime=float(findTime.match(line).group(1))
        if oldTime != -1:
            if oldTime != newTime:
                NCl.append(int(findClust.match(oldline).group(1)))
                Time.append(float(oldTime))
                oldTime=newTime
        else:
            oldTime=newTime
        oldline=line
NCl.append(int(findClust.match(oldline).group(1)))
Time.append(float(oldTime))


print("#%10s  %8s " % ("Time","Nclust"))
for n in range(len(Time)):
    print(" %10.2f  %8d " % (Time[n],NCl[n]+1))
    
