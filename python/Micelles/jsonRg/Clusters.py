'''
Created on Dec 23, 2017

@author: marchi
'''
import sys

import ujson as json
import numpy as np
from math import *
import collections as coll 

class Clusters(object):
    '''
    classdocs
    '''
    LARGEST_CLUSTER = 10
    @staticmethod
    def FilterOnSize(n):
        return lambda a, b: [h[1] for h in filter(lambda x: x[0] > n, [d for d in zip(a, b)])]

    def __init__(self, openfile=None, fileout=None, start_=None, end_=None, molecules=None):
        '''
        Constructor
        '''
        self.Rg = {'Pol': [], 'Tot': []}
        self.Rg_s = {'Pol': {}, 'Tot': {}}
        self.Wt = []
        self.Wt_s = {}
        self.ecc = {'Pol': [], 'Tot': []}
        self.ecc1 = {'Pol': [], 'Tot': []}
        self.ecc_s = {'Pol': {}, 'Tot': {}}
        self.ecc1_s = {'Pol': {}, 'Tot': {}}
        self.ax={'Pol': {}, 'Tot': {}}
        self.I={'Pol': {}, 'Tot': {}}
        self.aggr=coll.OrderedDict()
        self.aggrType=coll.OrderedDict()


        self.Start = 0
        self.End = 1e24
        self.start_ = start_
        self.end_ = end_
        if start_ is not None:
            self.Start = start_
        if end_ is not None:
            self.End = end_

        if fileout:
            self.fileout = open(fileout, 'w')
        else:
            self.fileout = sys.stdout

        if openfile:
            self.my_file = openfile
        else:
            print('\nNo stream given, abort \n')
            sys.exit(1)
        if not molecules:
            print('\nNo cluster composition, abort \n')
            sys.exit(1)
        else:
            self.molecules = molecules
            self.mols_first=molecules
        print(self.molecules)
        print(self.mols_first)
        sys.exit(1)

    def read(self):
        self.traj = json.load(self.my_file)
        self.List = [L for L in self.traj if float(L) >= self.Start and float(L) <= self.End]
        myStart = self.List[0]
        myEnd = self.List[len(self.List) - 1]
        if self.Start < float(myStart):
            self.Start = float(myStart)
        if self.End > float(myEnd):
            self.End = float(myEnd)

        if not len(self.List):
            print('Trajectory bounds wrong: Found %-s -- %-s, while given %10.3f -- %10.3f ' % (myStart, myEnd, self.Start, self.End))
        if self.start_ is None and self.end_ is None:
            print('   Iterating over the entire trajectory from step %-11.2f through %-11.2f ' % (self.Start, self.End))
        else:
            print('   Iterating from step %-11.2f through %-11.2f ' % (self.Start, self.End))

    def trajectory(self):
        def mymap(x):
            if x < 0:
                return 0
            else:
                return x

        Max=-10
        for time in self.List:
            clusters = self.traj[time]['cluster']
            N=0
            for cluster in range(len(clusters)):
                for Type in clusters[cluster]:
                    if Type in self.mols_first:
                        N+=1
            if N > Max:
                Max=N 
            
        for time in self.List:
            clusters = self.traj[time]['cluster']
            gyros = self.traj[time]['gyro']
            self.aggrType[time]=[0]*Max
            for Type in gyros:
                for gyro in range(len(gyros[Type])):
                    self.Rg[Type].append(gyros[Type][gyro]['Rg'])
                    if gyro not in self.Rg_s[Type]:
                        self.Rg_s[Type][gyro] = []
                        self.ax[Type][gyro] = []
                        self.I[Type][gyro] = []

                    self.Rg_s[Type][gyro].append(gyros[Type][gyro]['Rg'])

                    Inertia0 = gyros[Type][gyro]['I']
                    Inertia = [a  if a > 0 else 0 for a in Inertia0]
                    ax0 = gyros[Type][gyro]['ax']
                    ax = [a  if a > 0 else 0 for a in ax0]
                    self.ax[Type][gyro].append(ax)
                    self.I[Type][gyro].append(Inertia)

            self.aggr[time]=0
            for cluster in range(len(clusters)):
                if cluster not in self.Wt_s:
                    self.Wt_s[cluster] = {}

                nAtm = 0
                for Type in clusters[cluster]:
                    nAtm += clusters[cluster][Type][1]
                    if Type in self.mols_first:
                        Atmm=clusters[cluster][Type][0]
                        self.aggrType[time][cluster]=Atmm
                    if Type in self.molecules:
                        if Type not in self.Wt_s[cluster]:
                            self.Wt_s[cluster][Type] = []

                        self.Wt_s[cluster][Type].append(clusters[cluster][Type])
                self.Wt.append(nAtm)
                self.aggr[time]+=1 if nAtm > Clusters.LARGEST_CLUSTER else 0

        myLambda = Clusters.FilterOnSize(Clusters.LARGEST_CLUSTER)
        for Type in self.Rg:
            self.Rg[Type] = myLambda(self.Wt, self.Rg[Type])
        self.Wt = myLambda(self.Wt, self.Wt)

    def avg(self):
        sys.stdout.write(' Type      Rg-Zv         Rg-Zno        Rg-Avg         Ecc-In        Ecc-ax  \n')
        sys.stdout.write('-' * 76 + '\n')
        ecc1={'Pol': [], 'Tot': []}
        ecc={'Pol': [], 'Tot': []}
        for Type in self.ax:
            for n in self.ax[Type]:
                ax=np.average(self.ax[Type][n],0)
                if ax[1]-ax[0] > ax[2]-ax[1]:
                    axM=sum(ax[1:])*0.5
                    axmin=ax[0]
                else:
                    axmin=sum(ax[:2])*0.5
                    axM=ax[2]
                if axmin < 7:
                    ecc[Type].append(0)
                    ecc1[Type].append(0)
                    continue
                I=np.average(self.I[Type][n],0)
                Imin=min(I)
                Iavg=sum(I)/3.0
                
                ecc[Type].append(1 - Imin / Iavg)
                axavg = axmin * axmin / (axM * axM)
                ecc1[Type].append(sqrt(1-axavg))
        for Type in self.Rg:
            ww = list(map(lambda x: (4.0 * pi / 3.0) * x * x * x, self.Rg[Type]))
            myAvgww = np.average(self.Rg[Type], weights=ww)
            myAvgW = np.average(self.Rg[Type], weights=self.Wt)
            myEccW = np.average(ecc[Type],weights=[1 if a > 0 else 0 for a in ecc[Type]])
            myEcc1W = np.average(ecc1[Type],weights=[1 if a > 0 else 0 for a in ecc1[Type]])
            myAvg = np.average(self.Rg[Type])
            sys.stdout.write(' %-4s   ' % Type)
            sys.stdout.write('   %-10.3f    %-10.3f    %-10.3f   ' % (myAvgww, myAvgW, myAvg))
            sys.stdout.write('   %-10.3f    %-10.3f  ' % (myEccW, myEcc1W))
            sys.stdout.write(' \n')
        sys.stdout.write(' \n')
        sys.stdout.write(' \n')
        
        for Type in self.Rg:
            sys.stdout.write(' Type   Cluster No.      Rg-Avg          ax            ay            az          Ecc-In        Ecc-ax   \n')
            sys.stdout.write('-' * 104 + '\n')
            for n in self.Rg_s[Type]:
                Rg = self.Rg_s[Type][n]
                ax0=self.ax[Type][n]
                Rg_avg = np.average(Rg)
                if Rg_avg < 3.0:
                    continue
                ecc_avg = ecc[Type][n]
                ecc1_avg = ecc1[Type][n]
                ax=np.average(ax0,0)
                sys.stdout.write(' %-4s      %2d         ' % (Type, n))
                sys.stdout.write('   %-10.3f    %-10.3f    %-10.3f     %-10.3f    %-10.3f    %-10.3f  ' 
                                 % (Rg_avg, ax[0],ax[1],ax[2],ecc_avg, ecc1_avg))
                sys.stdout.write(' \n')
            sys.stdout.write(' \n')
        sys.stdout.write('                  Averaged Cluster Size        \n')
        sys.stdout.write('-' * 80 + '\n')
        for n in self.Wt_s:
            sys.stdout.write('%3d   ' % n)
            for Type in self.molecules:
                if Type in self.Wt_s[n]:
                    Wt = np.array(self.Wt_s[n][Type]).ravel()
                    Wt_avg0 = np.average(Wt[::3])
                    Wt_avg1 = np.average(Wt[1::3])
                    sys.stdout.write('%-s  %8.2f  %8.2f  ' % (Type, Wt_avg0, Wt_avg1))
                else:
                    sys.stdout.write('%-s  %8s  %8s  ' % (Type, (5 * ' ' + '-' + 2 * ' '), (5 * ' ' + '-' + 2 * ' ')))

            sys.stdout.write('\n')

    def histograms(self):
        histoRg = {}
#         histoEc = {}
#         histoEc1 = {}
        for Type in self.Rg:
            histoRg[Type] = np.histogram(self.Rg[Type], bins='auto', density=True)
#             histoEc[Type] = np.histogram(self.ecc[Type], bins='auto', density=True)
#             histoEc1[Type] = np.histogram(self.ecc1[Type], bins='auto', density=True)

        for Type in histoRg:
            y, x = histoRg[Type]
            self.fileout.write('# Rg %-s \n' % Type)
            for i in range(len(y)):
                self.fileout.write(' %10.2f  %11.5f \n' % (0.5 * (x[i] + x[i + 1]), y[i]))
            self.fileout.write('&\n')
#         for Type in histoEc:
#             y, x = histoEc[Type]
#             self.fileout.write('# Eccentricity %-s \n' % Type)
#             for i in range(len(y)):
#                 self.fileout.write(' %10.2f  %11.5f \n' % (0.5 * (x[i] + x[i + 1]), y[i]))
#             self.fileout.write('&\n')
#         for Type in histoEc1:
#             y, x = histoEc1[Type]
#             self.fileout.write('# Eccentricity  1 %-s \n' % Type)
#             for i in range(len(y)):
#                 self.fileout.write(' %10.2f  %11.5f \n' % (0.5 * (x[i] + x[i + 1]), y[i]))
#             self.fileout.write('&\n')
            
            
    def aggregation(self):
        lambda x,y:sum(y[:x])
        self.fileout.write('\n#  Time      No. Clusters  \n')
        for n in self.aggr:
            self.fileout.write('%-s  %4d ' % (n,self.aggr[n]))
            v=self.aggrType[n]
            values=[(lambda v,n:sum(v[:n+1]))(v,m) for m in range(len(v))]
            for m in range(len(v)):
                self.fileout.write('  %4d ' % v[m])
            self.fileout.write('\n' )
            
        

if __name__ == "__main__":
    filename = '/Users/marchi/AOT_RM10/tato3.json'
    fp = open(filename, 'r')
    mols = ['AOT', 'SOL', 'NA']
    clust = Clusters(openfile=fp, molecules=mols)
    clust.read()
    clust.trajectory()
    clust.avg()
    clust.histograms()
