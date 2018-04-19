'''
Created on Dec 23, 2017

@author: marchi
'''
import sys

import ujson as json
import numpy as np
from math import sqrt
from math import pi
import collections as coll
from . import PolyPep as pep

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
        self.Id=coll.OrderedDict()


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
            if 'Prot' in molecules:
                if type(molecules) is list:
                    molecules.remove('Prot')
                else:
                    molecules=[]
                for ami in pep.standard_aa_names:
                    molecules.append(ami)

            self.molecules = molecules
            self.mols_first=molecules[0]

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
        MolType=self.mols_first
        for time in self.List:
            clusters = self.traj[time]['cluster']
            gyros = self.traj[time]['gyro']
            for cluster in range(len(clusters)):
                if MolType not in clusters[cluster]['hsh']: continue
                myclust=clusters[cluster]['hsh'][MolType]
                if myclust not in self.Id:
                    self.Id[myclust]=[]
                self.Id[myclust].append(time)
                if myclust not in self.Wt_s:
                    self.Wt_s[myclust] = {}

            for Type in gyros:
                for gyro in range(len(gyros[Type])):
                    if MolType not in gyros[Type][gyro]['hsh']: continue
                    myclust=gyros[Type][gyro]['hsh'][MolType]
                    if myclust not in self.Rg_s[Type]:
                        self.Rg_s[Type][myclust] = []
                        self.ax[Type][myclust] = []
                        self.I[Type][myclust] = []

            ok=True
            for Type in gyros:
                for gyro in range(len(gyros[Type])):
                    if MolType not in gyros[Type][gyro]['hsh']: continue
                    self.Rg[Type].append(gyros[Type][gyro]['Rg'])
                    myclust=gyros[Type][gyro]['hsh'][MolType]
                    self.Rg_s[Type][myclust].append(gyros[Type][gyro]['Rg'])

                    Inertia0 = gyros[Type][gyro]['I']
                    Inertia = [a  if a > 0 else 0 for a in Inertia0]
                    ax0 = gyros[Type][gyro]['ax']
                    ax = [a  if a > 0 else 0 for a in ax0]
                    self.ax[Type][myclust].append(ax)
                    self.I[Type][myclust].append(Inertia)
                ok=False
            self.aggr[time]=0
            for cluster in range(len(clusters)):
                if MolType not in clusters[cluster]['hsh']: continue
                myclust=clusters[cluster]['hsh'][MolType]
                nAtm = 0
                for Type in clusters[cluster]:
                    if Type == 'hsh':
                        continue
                    nAtm += clusters[cluster][Type][1]
                    if Type in self.molecules:
                        if Type not in self.Wt_s[myclust]:
                            self.Wt_s[myclust][Type] = []
                        self.Wt_s[myclust][Type].append(clusters[cluster][Type][:2])
                self.Wt.append(nAtm)
                self.aggr[time]+=1 if nAtm > Clusters.LARGEST_CLUSTER else 0
        myLambda = Clusters.FilterOnSize(Clusters.LARGEST_CLUSTER)
        for Type in self.Rg:
            self.Rg[Type] = myLambda(self.Wt, self.Rg[Type])
        self.Wt = myLambda(self.Wt, self.Wt)

    def avg(self):
        sys.stdout.write(' Type      Rg-Zv         Rg-Zno        Rg-Avg         Ecc-In        Ecc-ax  \n')
        sys.stdout.write('-' * 76 + '\n')
        ecc1={'Pol': {}, 'Tot': {}}
        ecc={'Pol': {}, 'Tot': {}}
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
                    ecc[Type][n]=0
                    ecc1[Type][n]=0
                    continue
                I=np.average(self.I[Type][n],0)
                Imin=min(I)
                Iavg=sum(I)/3.0

                ecc[Type][n]=(1 - Imin / Iavg)
                axavg = axmin * axmin / (axM * axM)
                ecc1[Type][n]=(sqrt(1-axavg))

        for Type in self.Rg:
            myEcc=[ecc[Type][n] for n in ecc[Type]]
            myEcc1=[ecc1[Type][n] for n in ecc1[Type]]
            ww = list(map(lambda x: (4.0 * pi / 3.0) * x * x * x, self.Rg[Type]))
            myAvgww = np.average(self.Rg[Type], weights=ww)
            myAvgW = np.average(self.Rg[Type], weights=self.Wt)
            myEccW = np.average(myEcc,weights=[1 if a > 0 else 0 for a in myEcc])
            myEcc1W = np.average(myEcc1,weights=[1 if a > 0 else 0 for a in myEcc1])
            myAvg = np.average(self.Rg[Type])
            sys.stdout.write(' %-4s   ' % Type)
            sys.stdout.write('   %-10.3f    %-10.3f    %-10.3f   ' % (myAvgww, myAvgW, myAvg))
            sys.stdout.write('   %-10.3f    %-10.3f  ' % (myEccW, myEcc1W))
            sys.stdout.write(' \n')
        sys.stdout.write(' \n')
        sys.stdout.write(' \n')

        for Type in self.Rg:
            sys.stdout.write(' Type   Cluster No.      Rg-Avg          ax            ay            az          Ecc-In        Ecc-ax      Occurence\n')
            sys.stdout.write('-' * 118 + '\n')
            N=0
            for n in self.Rg_s[Type]:
                Rg = self.Rg_s[Type][n]
                ax0=self.ax[Type][n]
                Rg_avg = np.average(Rg)
                if Rg_avg < 3.0:
                    continue
                ecc_avg = ecc[Type][n]
                ecc1_avg = ecc1[Type][n]
                ax=np.average(ax0,0)
                sys.stdout.write(' %-4s      %2d         ' % (Type, N))
                sys.stdout.write('   %-10.3f    %-10.3f    %-10.3f     %-10.3f    %-10.3f    %-10.3f   %d'
                                 % (Rg_avg, ax[0],ax[1],ax[2],ecc_avg, ecc1_avg,len(Rg)))
                sys.stdout.write(' \n')
                N+=1
            sys.stdout.write(' \n')
        sys.stdout.write('                  Averaged Cluster Size        \n')
        sys.stdout.write('-' * 80 + '\n')
        N=0
        for n in self.Wt_s:
            sys.stdout.write('%3d   ' % N)
            for Type in self.molecules:
                if Type in self.Wt_s[n]:
                    Wt = np.array(self.Wt_s[n][Type]).ravel()
                    Wt_avg0 = np.average(Wt[::2])
                    sys.stdout.write('%-s  %8.2f  ' % (Type, Wt_avg0))
            sys.stdout.write('   %9d  ' % len(Wt[::2]))
            sys.stdout.write('\n')
            N+=1
    def histograms(self):
        histoRg = {}
        for Type in self.Rg:
            histoRg[Type] = np.histogram(self.Rg[Type], bins='auto', density=True)

        for Type in histoRg:
            y, x = histoRg[Type]
            self.fileout.write('# Rg %-s \n' % Type)
            for i in range(len(y)):
                self.fileout.write(' %10.2f  %11.5f \n' % (0.5 * (x[i] + x[i + 1]), y[i]))
            self.fileout.write('&\n')

    def aggregation(self):
        if self.fileout is sys.stdout:
            print('\n\n         Need output file to plot aggregation       \n')
            print('!'*90)
            sys.exit(1)
        uscores=['_Ag','_Rg','_E0','_E1','_Cl']
        fn=self.fileout.name
        filenames=[]
        if fn.find('.') >0:
            for uscore in uscores:
                filenames+=[fn[:fn.find('.')]+uscore+fn[fn.find('.'):]]
        else:
            for uscore in uscores:
                filenames+=[fn+uscore]
        fps=[open(fname,'w') for fname in filenames]
        for fp in fps:
            if 'Cl' in fp.name:
                continue
            M=0
            for id0 in self.Wt_s:
                fp.write('@    s%-d legend "Id No. %-s"\n' % (M,M))
                fp.write('@    s%-d symbol 1 \n' % (M))
                fp.write('@    s%-d symbol size 0.35 \n' % (M))
                fp.write('@    s%-d symbol fill pattern 1\n' % (M))
                fp.write('@    s%-d line type 0 \n' % (M))
                M+=1
        for id0 in self.Wt_s:
            Type=next(iter(self.Wt_s[id0]))
            n=0
            for time in self.Id[id0]:
                ax=self.ax['Tot'][id0][n]
                Rg=self.Rg_s['Tot'][id0][n]
                I=self.I['Tot'][id0][n]
                naggr=self.Wt_s[id0][Type][n][0]
                if ax[1]-ax[0] > ax[2]-ax[1]:
                    axM=sum(ax[1:])*0.5
                    axmin=ax[0]
                else:
                    axmin=sum(ax[:2])*0.5
                    axM=ax[2]
                if axmin < 7:
                    ecc=0
                    ecc1=0
                else:
                    Imin=min(I)
                    Iavg=sum(I)/3.0
                    ecc=(1 - Imin / Iavg)
                    axavg = axmin * axmin / (axM * axM)
                    ecc1=(sqrt(1-axavg))
                if Rg > 5:
                    fps[0].write('%-s  %4d\n'% (time,naggr))
                    fps[1].write('%-s  %10.3f\n'% (time,Rg))
                    fps[2].write('%-s  %10.3f\n'% (time,ecc))
                    fps[3].write('%-s  %10.3f\n'% (time,ecc1))
                n+=1
            [fps[n].write('&\n') for n in range(4)]

        for n in self.aggr:
            fps[4].write('%-s  %4d \n' % (n,self.aggr[n]))
        [fp.close() for fp in fps]



if __name__ == "__main__":
    filename = '/Users/marchi/AOT_RM10/mytest.json'
    fp = open(filename, 'r')
    mols = ['AOT', 'SOL', 'NA','Prot']
    clust = Clusters(openfile=fp, molecules=mols)
    clust.read()
    clust.trajectory()
    clust.avg()
    clust.histograms()
