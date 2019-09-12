#!/Users/marchi/anaconda3/bin/python3
import re
import sys
import argparse
import math
exclusion=['SOD','CLA','CAL','MG','ZN']

parser = argparse.ArgumentParser(prog='createNonnonded',description='Compute nonbonded weighted interactions')
parser.add_argument('--input','-i', action = 'store', type = str, help = 'Input pdb file', required=True)
parser.add_argument('--factor','-f', action = 'store', type = float, help = 'Fudge factor', default=1.1)
parser.add_argument('--water','-w', action = 'store', type = str, help = 'Water oxygen', default='OWT3')

args = parser.parse_args()
try:
    fp=open(args.input,'r')
except FileNotFoundError:
    print(args.input+": File not found!")
    sys.exit(2)
factw=args.factor;
typew=args.water;
lines=fp.readlines()
findLabel=re.compile(r'^ *\[  *(\w+)  *\].*')
findComment=re.compile(r'(.*);.*')
potential={}
label='null'
ok=True
for line in lines:
    if len(line.split()) == 0: continue
    ncmt=line.find(';')
    if ncmt >= 0:
        line=line[0:ncmt]
        if not len(line):
            continue
    lbl=findLabel.match(line)
    if lbl:
        label=lbl.group(1)
        potential[label]=[]
        continue
    if re.search(r'.*HEAVY_H.*',line):
        ok=False
    if re.search(r'#else.*',line):
        ok=True
    if re.search(r'#endi.*',line):
        ok=True
    if ok:
        if re.search(r'^#.*',line):
            continue
        potential[label].append(line)

key='atomtypes'
potent=potential['atomtypes']
pots=[]
for line in potent:
    tokens=line.split()
    sigma=float(tokens[5])
    eps=float(tokens[6])
    pot0=[tokens[0],sigma,eps]
    pots.append(pot0)

print('[ pairtypes ]\n; i	j	func	sigma1-4	epsilon1-4 ; THESE ARE 1-4 INTERACTIONS')
ppots={}
for i in range(len(pots)):
    t_i=pots[i][0]
    sigma_i=pots[i][1]
    eps_i=pots[i][2]
    ppots[t_i]={}
    for j in range(i,len(pots)):
        t_j=pots[j][0]
        sigma_j=pots[j][1]
        eps_j=pots[j][2]
        sigma=0.5*(sigma_i+sigma_j)
        eps=math.sqrt(eps_i*eps_j)
        ppots[t_i][t_j]=[sigma,eps]
        if t_i.strip() == typew.strip() and t_j.strip() == typew.strip():
            continue
        if t_i.strip() in exclusion or t_j.strip() in exclusion:
            continue
        if t_i.strip() == typew.strip() or t_j.strip() == typew.strip():
            ppots[t_i][t_j]=[sigma,eps*factw]
        
        

key='pairtypes'
potent=potential[key]
for line in potent:
    tokens=line.split()
    t_1=tokens[0]
    t_4=tokens[1]
    sigma=float(tokens[3])
    eps=float(tokens[4])
    if t_1 in ppots:
        if t_4 in ppots[t_1]:
            ppots[t_1][t_4]=[sigma,eps]
    if t_4 in ppots:
        if t_1 in ppots[t_4]:
            ppots[t_4][t_1]=[sigma,eps]
        

for t_1 in ppots:
    for t_4 in ppots[t_1]:
        sigma=ppots[t_1][t_4][0]
        eps=ppots[t_1][t_4][1]
        print("%-5s %-5s  1  %16.12f  %16.12f " % (t_1, t_4, sigma, eps))

print('[ nonbond_params ]')
print('  ; i    j func       V(sigma)        W(epsilon)')
for i in range(len(pots)):
    t_i=pots[i][0]
    for j in range(i,len(pots)):
        t_j=pots[j][0]
        sigma=ppots[t_i][t_j][0]
        eps=ppots[t_i][t_j][1]
        print("%-5s %-5s  1  %16.12f  %16.12f " % (t_i, t_j, sigma, eps))

