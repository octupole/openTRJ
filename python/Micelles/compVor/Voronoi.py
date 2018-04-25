'''
Created on Dec 15, 2017

@author: marchi
'''
import sys
import ujson as json
import Micelles.compVor.jsonIterator as jsonIterator
import numpy as np
import collections as coll


_Res='Res'
_List='List'
_Types='Types'
_Cofact='Cofact'
_Water='Water'
_Prot='Prot'
_Ions='Ions'
_Surf='Surf'
_Other='Other'
_OrgSol='OrgSol'
_Dna='Dna'
_DetPol='DetPol'
_DetOil='DetOil'
_Area='Area'

class Voronoi(object):
    '''
    classdocs
    '''
    Labels=['Vol','Shell','AClust','Pick']
    def __init__(self, openfile=None,fileout=None,start=None,end=None):
        '''
        Constructor
        '''
        funcs=[self.avgVol,self.avgShell,self.avgAClust,self.pickRes]
        self.whatToDo=dict(zip(Voronoi.Labels,funcs))
        self.start=None
        self.end=None
        if start:
            self.start=float(start)
        if end:
            self.end=float(end)

        if fileout:
            self.fileout=open(fileout,'w')
        else:
            self.fileout=sys.stdout

        if openfile:
            self.my_file=openfile
        else:
            print('\nNo stream given, abort \n')
            sys.exit(1)
    @classmethod
    def testLabels(cls,label):
        if label in cls.Labels:
            return True
        return False

    def what(self,label):
        if label in self.whatToDo:
            return self.whatToDo[label]
        return None

    def read(self):
        self.root=json.load(self.my_file)
        print('Read in')
        self.Atom=self.root['Res']['List']
        self.myIter=jsonIterator.JSONIterator(self.root,self.start,self.end)

    def avgVol(self):
        root=self.root
        myIter=self.myIter
        Atom=self.Atom

        nType={}
        volTot={}
        timeC=myIter.next()
        while timeC != None:
            vols=root[timeC]['Vol']
            n=0
            for vol in vols:
                typeA=Atom[n]
                if typeA not in volTot:
                    volTot[typeA]=[]
                    nType[typeA]=0
                volTot[typeA].append(vol)
                nType[typeA]+=1
                n+=1
            timeC=myIter.next()
        for key in volTot:
            print(' %s  %10.3f %10.5f '% (key,np.average(volTot[key]), np.std(volTot[key])))

    def avgShell(self):
        root=self.root
        myIter=self.myIter
        timeC=myIter.next()
        shell={}
        while timeC != None:
            vols=root[timeC]['Shell']
            if not isinstance(vols, dict):
                print("\n--------- Cannot find Shell in json input, stop here  -------------\n")
                sys.exit()
                
                
            for vol in vols:
                if vol not in shell:
                    shell[vol]={'No': [],'Vol': []}
                shell[vol]['No'].append(vols[vol]['No'])
                shell[vol]['Vol'].append(vols[vol]['Vol'])
            timeC=myIter.next()

        sys.stdout.write('Level    No Wat     std          Vol/mol     std\n')
        for keys in shell:
            sys.stdout.write('  %s '% keys)
            for key in shell[keys]:
                sys.stdout.write(' %10.2f %10.4f '%
                                 (np.average(shell[keys][key]),np.std(shell[keys][key])))
            sys.stdout.write('\n')

    def writePick(self,histS):
        n=0
        self.fileout.write('  No  Residue  ')
        myhistS=next(iter(histS.values()))
        for label in myhistS:
            self.fileout.write('      %-7s     '%label)
        self.fileout.write('\n')
        for resn in self.myPick:
            res=self.root[_Res][_List][resn]
            hist=histS[resn]
            self.fileout.write('%3d     %-3s    '% (resn,res))
            for label in hist:
                avg=np.average(hist[label])
                std=np.std(hist[label])
#                print(len(hist[label]))
                self.fileout.write(' %7.3f (%6.4f) '% (avg,std))
            self.fileout.write('\n')
            n+=1

    def doTypes(self,myPick):
        self.myPick=myPick
        root=self.root
        myIter=self.myIter
        timeC=myIter.next()
        histS=coll.OrderedDict()
        totAreaS=coll.OrderedDict()
        while timeC != None:
            for resn in myPick:
                res=self.root[_Res][_List][resn]
                if res not in histS:
                    histS[res]={}
                    totAreaS[res]=[]
                hist=histS[res]
                totArea=totAreaS[res]
                area=root[timeC][_Area][resn]
                tot=0.0
                for _type in area:
                    tot+=area[_type]

                for _type in area:
                    if _type not in hist:
                        hist[_type]=[]
                    hist[_type].append(area[_type]/tot)
                totArea.append(tot)
            timeC=myIter.next()
        self.fileout.write('   Residue ')
        for label in hist:
            self.fileout.write('      %-7s     '%label)
        self.fileout.write('\n')
        for res in histS:
            self.fileout.write('    %-4s   '%res)
            for label in hist:
                avg=np.average(histS[res][label])
                std=np.std(histS[res][label])
                self.fileout.write(' %7.3f (%6.4f) '% (avg,std))
            self.fileout.write('\n')

    def doLists(self,myPick):
        self.myPick=myPick
        root=self.root
        myIter=self.myIter
        timeC=myIter.next()
        histS=coll.OrderedDict()
        totAreaS=coll.OrderedDict()
        while timeC != None:
            for resn in myPick:
                if resn not in histS:
                    histS[resn]={}
                    totAreaS[resn]=[]
                hist=histS[resn]
                totArea=totAreaS[resn]
                area=root[timeC][_Area][resn]
                tot=0.0
                for _type in area:
                    tot+=area[_type]

                for _type in area:
                    if _type not in hist:
                        hist[_type]=[]
                    hist[_type].append(area[_type]/tot)
                totArea.append(tot)
            timeC=myIter.next()
        self.histS=histS

    def pickRes(self,myPick):

        
        root=self.root
        if isinstance(myPick,list):
            if isinstance(myPick[0],int):
                self.doLists(myPick)
                self.writePick(self.histS)
            elif isinstance(myPick[0],str):
                res0=myPick
                print(res0)
                res=[]
                res_n=[]
                for res1 in res0:
                    ok=False
                    if res1 in root[_Res]:
                        ok=True
                        res.append(res1)
                    if res1+'_O' in root[_Res]:
                        ok=True
                        res.append(res1+'_O')
                    if  res1+'_P' in root[_Res]:
                        ok=True
                        res.append(res1+'_P')
                    if not ok:
                        res_n.append(res1)

                if len(res_n):
                    print('Wrong residue entered. Some residues do not exist:')
                    for res1 in res_n:
                        print('Unfound %-s'% res1)
                    sys.exit(1)
                indexPicks=[i for i in range(len(root[_Res][_List])) if root[_Res][_List][i] in res]
                self.doTypes(indexPicks)
        else:
            print('pickRes: Should enter an int, a str or a list of in or str instead entered a %-s ' % type(myPick))
            sys.exit(-1)
    def avgAClust(self):
        root=self.root
        myIter=self.myIter
        timeC=myIter.next()
        area={}

        while timeC != None:
            clusters=root[timeC]['AClust']
            vClusters=root[timeC]['Clust']
            for cluster in range(len(clusters)):
                if clusters[cluster]['NoAtm'] < 200:
                    continue
                if cluster not in area:
                    area[cluster]={}
                for elem in clusters[cluster]:
                    if elem not in area[cluster]:
                        area[cluster][elem]=[]
                    area[cluster][elem].append(clusters[cluster][elem])
                if 'vol' not in area[cluster]:
                    area[cluster]['vol']=[]
                area[cluster]['vol'].append(vClusters[cluster])
            timeC=myIter.next()

        for cluster in area:
            self.fileout.write('\n         Cluster No.    %s \n'% cluster)
            self.fileout.write('%-10s  %10s %10s \n'% ('  Type ',' Surface ','  std '))
            tmp=[key for key in area[cluster] if key not in ['Tot','NoAtm','vol']]+['Tot','NoAtm','vol']
            for key in tmp:
                avg=np.average(area[cluster][key])
                std=np.std(area[cluster][key])
                self.fileout.write('%-10s  %10.2f %10.4f \n'% (key,avg,std))




if __name__ == "__main__":
    import time
    filename='/Users/marchi/AOTs/2AML/t600-1200-det.json'
    fp=open(filename,'r')
    myVor=Voronoi(fp)
    start_time = time.time()
    myVor.read()
    print("--- %s seconds ---" % (time.time() - start_time))
    start_time = time.time()
    myVor.pickRes(list(range(40)))
    print("--- %s seconds ---" % (time.time() - start_time))
    del myVor

