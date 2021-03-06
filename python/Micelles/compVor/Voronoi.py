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
_Vol='Vol'

class Voronoi(object):
    '''
    classdocs
    '''
    Labels=['Vol','Shell','AClust','Pick','TotVol']
    def __init__(self, openfile=None,fileout=None,start=None,end=None,type='seq'):
        '''
        Constructor
        '''
        funcs=[self.avgVol,self.avgShell,self.avgAClust,self.pickRes,self.TotVol]

        self.whatToDo=dict(zip(Voronoi.Labels,funcs))

        self.start=None
        self.end=None
        self.type=type
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
        self.myIter=jsonIterator.JSONIterator(self.my_file,self.start,self.end,self.type)
#        self.Atom=json.loads(json.dumps(self.myIter['Res']['List']))
        self.AtomTypes=json.loads(json.dumps(self.myIter['Res']))
        self.Atom=self.AtomTypes['List']

    def TotVol(self):
        myIter=self.myIter
        Atom=self.Atom

        volTot={}
        natoTot={}
        for type in self.AtomTypes:
            if type != 'List':
                typeB=self.AtomTypes[type][1]
                volTot[typeB]=[]
                natoTot[typeB]=0
        types=[type for type in volTot]
        timeC=next(myIter)
        firstRound=True
        while timeC != None:
            vols=myIter[timeC]['Vol']
            volTypes={}
            for type in types:
                volTypes[type]=[]
            n=0
            for vol in vols:
                typeA=Atom[n]
                typeB=self.AtomTypes[typeA][1]
                volTypes[typeB].append(vol)
                n+=1
            if firstRound:
                firstRound=False
                n=0
                for vol in vols:
                    typeA=Atom[n]
                    typeB=self.AtomTypes[typeA][1]
                    natoTot[typeB]+=1
                    n+=1

            for type in types:
                Tot=np.sum(volTypes[type])
                volTot[type].append(Tot)
            
            timeC=next(myIter)
        for key in volTot:
            print(' %-6s  %10.3f %10.5f   natom = %7d'% (key,np.average(volTot[key]), np.std(volTot[key]),natoTot[key]))

    def avgVol(self):
        myIter=self.myIter
        Atom=self.Atom

        nType={}
        volTot={}
        timeC=next(myIter)
        while timeC != None:
            vols=myIter[timeC]['Vol']
            n=0
            for vol in vols:
                typeA=Atom[n]
                if typeA not in volTot:
                    volTot[typeA]=[]
                    nType[typeA]=0
                volTot[typeA].append(vol)
                nType[typeA]+=1
                n+=1
            timeC=next(myIter)
        for key in volTot:
            print(' %s  %10.3f %10.5f '% (key,np.average(volTot[key]), np.std(volTot[key])))


    def avgShell(self):
        myIter=self.myIter
        timeC=next(myIter)
        shell={}
        while timeC != None:
            vols=myIter[timeC]['Shell']
            if not isinstance(vols, dict):
                print("\n--------- Cannot find Shell in json input, stop here  -------------\n")
                sys.exit()
                
                
            for vol in vols:
                if vol not in shell:
                    shell[vol]={'No': [],'Vol': []}
                shell[vol]['No'].append(vols[vol]['No'])
                shell[vol]['Vol'].append(vols[vol]['Vol'])
            timeC=next(myIter)

        sys.stdout.write('Level    No Wat     std          Vol/mol     std\n')
        for keys in shell:
            sys.stdout.write('  %s '% keys)
            for key in shell[keys]:
                sys.stdout.write(' %10.2f %10.4f '%
                                 (np.average(shell[keys][key]),np.std(shell[keys][key])))
            sys.stdout.write('\n')

    def writePick(self):
        histS=self.histS
        totAreaS=self.totAreaS
        n=0
        self.fileout.write('  No  Residue  ')
        myhistS=next(iter(histS.values()))
        self.fileout.write('      %-7s     '%'Vol')
        for label in myhistS:
            self.fileout.write('      %-7s Area'%label)
        self.fileout.write('       %-10s     '%'Total Area')            
        self.fileout.write('\n')
        for resn in self.myPick:
            res=self.Atom[resn]
            hist=histS[resn]
            totArea=totAreaS[resn]

            self.fileout.write('%3d     %-3s    '% (resn,res))
            totVol=self.totVolS[resn]
            avg=np.average(totVol)
            std=np.std(totVol)
            self.fileout.write('   %7.3f (%6.4f) '% (avg,std))

            for label in hist:
                avg=np.average(hist[label])
                std=np.std(hist[label])
                self.fileout.write(' %7.3f (%6.4f) '% (avg,std))

                
            avg=np.average(totArea)
            std=np.std(totArea)
            self.fileout.write('   %7.3f (%6.4f) '% (avg,std))
            self.fileout.write('\n')
            n+=1

    def doTypes(self,myPick):
        self.myPick=myPick
        myIter=self.myIter
        timeC=next(myIter)
        histS=coll.OrderedDict()
        totAreaS=coll.OrderedDict()
        while timeC != None:
            for resn in myPick:
                res=self.Atom[resn]
                if res not in histS:
                    histS[res]={}
                    totAreaS[res]=[]
                hist=histS[res]
                totArea=totAreaS[res]
                area=self.myIter[timeC][_Area][resn]
                tot=0.0
                for _type in area:
                    tot+=area[_type]

                for _type in area:
                    if _type not in hist:
                        hist[_type]=[]
                    hist[_type].append(area[_type]/tot)
                totArea.append(tot)
            timeC=next(myIter)
        self.fileout.write('   Residue ')
        for label in hist:
            self.fileout.write('      %-7s     '%label)
        
        self.fileout.write('       %-10s     '%'Total Area')
        self.fileout.write('\n')
        for res in histS:
            self.fileout.write('    %-4s   '%res)
            for label in hist:
                avg=np.average(histS[res][label])
                std=np.std(histS[res][label])
                self.fileout.write(' %7.3f (%6.4f) '% (avg,std))
            avg=np.average(totAreaS[res])
            std=np.std(totAreaS[res])
            self.fileout.write('  %10.3f (%6.6f)' % (avg,std))
            self.fileout.write('\n')
        

    def doLists(self,myPick):
        self.myPick=myPick
        myIter=self.myIter
        timeC=next(myIter)
        histS=coll.OrderedDict()
        totAreaS=coll.OrderedDict()
        totVolS=coll.OrderedDict()
        while timeC != None:
            for resn in myPick:
                if resn not in histS:
                    histS[resn]={}
                    totAreaS[resn]=[]
                    totVolS[resn]=[]
                hist=histS[resn]
                totArea=totAreaS[resn]
                totVol=totVolS[resn]
                area=myIter[timeC][_Area][resn]
                totVol.append(myIter[timeC][_Vol][resn])
                tot=0.0
                for _type in area:
                    tot+=area[_type]

                for _type in area:
                    if _type not in hist:
                        hist[_type]=[]
                    hist[_type].append(area[_type]/tot)
                totArea.append(tot)
            timeC=next(myIter)
        self.histS=histS
        self.totAreaS=totAreaS
        self.totVolS=totVolS

    def pickRes(self,myPick):
        myIter=self.myIter
        if isinstance(myPick,list):
            if isinstance(myPick[0],int):
                self.doLists(myPick)
                self.writePick()
            elif isinstance(myPick[0],str):
                res0=myPick
                res=[]
                res_n=[]
                for res1 in res0:
                    ok=False
                    if res1 in myIter[_Res]:
                        ok=True
                        res.append(res1)
                    if res1+'_O' in myIter[_Res]:
                        ok=True
                        res.append(res1+'_O')
                    if  res1+'_P' in myIter[_Res]:
                        ok=True
                        res.append(res1+'_P')
                    if not ok:
                        res_n.append(res1)

                if len(res_n):
                    print('Wrong residue entered. Some residues do not exist:')
                    for res1 in res_n:
                        print('Unfound %-s'% res1)
                    sys.exit(1)
                indexPicks=[i for i in range(len(myIter[_Res][_List])) if myIter[_Res][_List][i] in res]
                self.doTypes(indexPicks)
        else:
            print('pickRes: Should enter an int, a str or a list of in or str instead entered a %-s ' % type(myPick))
            sys.exit(-1)

    def avgAClust(self):
        myIter=self.myIter
        timeC=next(myIter)
        area={}

        while timeC != None:
            clusters=myIter[timeC]['AClust']
            vClusters=myIter[timeC]['Clust']
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
            timeC=next(myIter)

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

