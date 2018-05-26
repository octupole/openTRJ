import sys

def replaceI(filein,fileout):
    fp=open(filein,'r')
    fpo=open(fileout,'w')
    N=13
    for line in fp:
        newLine=line
        if line.find('C1') > -1 or line.find('C6') > -1 or line.find('C7') > -1 or line.find('C8') > -1 or line.find('C5') > -1:
            newLine=line[:N]+'Iz'+line[N+2:]
        elif line.find('C2') > -1:
            newLine=line[:N]+'C '+line[N+2:]
        elif line.find('C3') > -1:
            newLine=line[:N]+'Iy'+line[N+2:]
        elif line.find('C4') > -1:
            newLine=line[:N]+'Ix'+line[N+2:]
        fpo.write(newLine)

def replaceU(filein,fileout):
    fp=open(filein,'r')
    fpo=open(fileout,'w')
    N=13
    for line in fp:
        newLine=line
        if line.find('C1') > -1 or line.find('C6') > -1 or line.find('C7') > -1 or line.find('C8') > -1 or line.find('C5') > -1:
            newLine=line[:N]+'Uz'+line[N+2:]
        elif line.find('C2') > -1:
            newLine=line[:N]+'C '+line[N+2:]
        elif line.find('C3') > -1:
            newLine=line[:N]+'Uy'+line[N+2:]
        elif line.find('C4') > -1:
            newLine=line[:N]+'Ux'+line[N+2:]
        fpo.write(newLine)        

if __name__ == '__main__':
    fileout='FILE.pdb'
    what='I'
    if len(sys.argv) < 2:
        print('at least an argument needed for '+sys.argv[0])
        sys.exit(1)
    elif len(sys.argv) > 4:
        print('Not more than three arguments are neede by '+sys.argv[0])
        sys.exit(1)
    elif len(sys.argv) == 3:
        fileout=sys.argv[2]
    elif len(sys.argv) == 4:
        fileout=sys.argv[2]
        what=sys.argv[3]

    filein=sys.argv[1]
    if what == 'I':
        replaceI(filein,fileout)
    else:
        replaceU(filein,fileout)
        