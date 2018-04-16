'''
Created on Dec 14, 2017

@author: marchi
'''
import re
import sys
import Micelles.Utilities.Averages as av

class Clusters(object):
    '''
    classdocs
    '''
    reClust=re.compile(' *Cluster')
    reTime=re.compile('  *.*time =(..*) *Rg')
    myClt={}
    start=0
    end=1e20

    def __init__(self,openfile=None,fileout=None,start=None,end=None,molecules=None):
        '''
        Constructor
        '''
        if start:
            self.start=float(start)
        if end:
            self.end=float(end);
            
        if fileout:
            self.fileout=open(fileout,'w')
        else:
            self.fileout=sys.stdout

        if openfile:
            self.my_file=openfile
        else:
            print('\nNo stream given, abort \n')
            sys.exit(1)
        if not molecules:
            print('\nNo cluster composition, abort \n')
            sys.exit(1)
        
        self.reMols=[]
        for key in molecules:
            reKey=re.compile(' *Cluster .*'+key+'\[(\d\d*)\]')
            self.reMols.append([key,reKey])
        
    def read(self):
        n=0
        for line in self.my_file:
            Time=re.match(self.reTime, line)
            if Time:
                time=re.match(self.reTime, line).groups()[0].strip()
            if re.match(self.reClust,line):
                for ree in self.reMols:
                    Match=re.match(ree[1],line)
                    if Match:
                        if time not in self.myClt:
                            self.myClt[time]={}
                            print(time)
                        num=float(Match.groups()[0].strip())
                    if float(time) > self.start and float(time) < self.end:
                        if ree[0] not in self.myClt[time]:
                            self.myClt[time][ree[0]]=[]
                        self.myClt[time][ree[0]].append(num)
    
    def avg(self):
        if len(self.myClt) == 0:
            print("Must accumulate cluster sizes first !")
            sys.exit(1)
            
           
        avgCL={}
        for ree in self.reMols:
            avgCL[ree[0]]=av.Averages()
        
#         for Cl in self.myClt:
#             self.fileout.write(' %s  '% Cl)
#             for key,item in self.myClt[Cl].items():
#                 self.fileout.write(' %3d  ' % len(item))
#             self.fileout.write(' \n')
        for Cl in self.myClt:
            gig=[]
            avi=0
            N=0
            for key,item in self.myClt[Cl].items():
                gig.append(item)
                avi+=len(item)
                N+=1
            avi/=N
            avf=0
            N=0
            for it in range(len(gig[0])):
                avf+=gig[1][it]/gig[0][it]
                N+=1
            avf/=N
            self.fileout.write(' %s  %2d %6.2f \n' % (Cl,avi,avf))
            

if __name__ == "__main__":
    filename='/Users/marchi/AOT_RM10/data.gyr'
    fp=open(filename,'r')
    mols=['AOT','SOL']
    clust=Clusters(openfile=fp,molecules=mols)
    clust.read()
    clust.avg()
    