'''
Created on Dec 12, 2017

@author: marchi
'''
import re
from math import *
import sys
import Micelles.Utilities.Averages as av


class Rg(object):
    reMic=re.compile(' *ClusMIC')
    rePol=re.compile(' *ClusPOL')
    reRg=re.compile('  *.*Rg =(..*) a =')
    reTime=re.compile('  *.*time =(..*) *Rg')
    myRgs={}
    myRgsP={}
    myRH={}
    myRHP={}
    myClt={}
    start=0
    end=1e16
    reClust=re.compile(' *Cluster .*Size (\d\d*) ')   
    
    def __init__(self,openfile=None,fileout=None,start=None,end=None):
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

    def read(self):
        n=0
        for line in self.my_file:
            if re.match(self.reMic,line) or re.match(self.rePol,line) :
                myline=line.split()
                a=float(myline[11])
                b=float(myline[14])
                c=float(myline[17])
                RH=pow(a*b*c,1.0/3.0)
                if isnan(RH):
                    RH=0
                Rg=re.match(self.reRg, line).groups()[0].strip()
                time=re.match(self.reTime, line).groups()[0].strip()
                if float(time) > self.end:
                    break
                if re.match(self.reMic,line) :
                    if float(time) > self.start and float(time) < self.end:
                        if time not in self.myRgs:
                            self.myRgs[time]=[]
                            self.myRH[time]=[]
                        if float(Rg) > 2.0:
                            self.myRgs[time].append(float(Rg))
                            self.myRH[time].append(float(RH))
                        if len(self.myRgs[time]) == 1:
                            if not n%100:
                                print(" %10.1f %10.2f %10.2f " % (float(time), float(Rg), RH))
                            n+=1
                else:
                    if float(time) > self.start and float(time) < self.end:
                        if time not in self.myRgsP:
                            self.myRgsP[time]=[]
                            self.myRHP[time]=[]

                        self.myRgsP[time].append(float(Rg))
                        self.myRHP[time].append(float(RH))
                        if len(self.myRgsP[time]) == 1:
                            if not n%100:
                                print(" %10.1f %10.2f %10.2f" % (float(time), float(Rg),RH))
                            n+=1
            clustMatch=re.match(self.reClust,line)
            if clustMatch:
                if float(time) > self.start and float(time) < self.end:
                    if time not in self.myClt:
                        self.myClt[time]=[]
                    num=float(clustMatch.groups()[0])
                    self.myClt[time].append(num)
                    


    def avg(self):
        if len(self.myRgs) == 0:
            print("Must accumulate Rg's first !")
            sys.exit(1)
            
        def Volume(x): return x*x*x
        avgRg0=av.Averages();
        avgRgV=av.Averages(Volume);
        avgRg0P=av.Averages();
        avgRgVP=av.Averages(Volume);
        avgRh0=av.Averages();
        avgRhV=av.Averages(Volume);
        avgRh0P=av.Averages();
        avgRhVP=av.Averages(Volume);
    
        for Rgs in self.myRgs:
            lavgRg=av.Averages();
            lavgRH=av.Averages();
            for n in range(len(self.myRgs[Rgs])):
                Rg=self.myRgs[Rgs][n]
                W=self.myClt[Rgs][n]
                
                RgP=self.myRgsP[Rgs][n]
                RH=self.myRH[Rgs][n]
                RHP=self.myRHP[Rgs][n]
                lavgRg.append(Rg,W)
                lavgRH.append(RH,W)
                avgRg0.append(Rg)
                avgRgV.append(Rg,W)
                avgRg0P.append(RgP)
                avgRgVP.append(RgP,W)
                avgRh0.append(RH)
                avgRhV.append(RH,W)
                avgRh0P.append(RHP)
                avgRhVP.append(RHP,W)
            self.fileout.write(' %10.1f %10.3f %10.3f \n' % (float(Rgs), lavgRg(),lavgRH()))
        print(' avgRg  = %10.3f avgRgV  = %10.3f ' % (avgRg0(),avgRgV()))
        print(' avgRgP = %10.3f avgRgVP = %10.3f ' % (avgRg0P(),avgRgVP()))
        print(' avgRh  = %10.3f avgRhV  = %10.3f ' % (avgRh0(),avgRhV()))
        print(' avgRhP = %10.3f avgRhVP = %10.3f ' % (avgRh0P(),avgRhVP()))

        