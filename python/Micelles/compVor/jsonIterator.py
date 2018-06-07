'''
Created on Dec 17, 2017

@author: marchi
'''
import sys 
import io
import json
import subprocess
class JSONIterator(object):
    '''
    classdocs
    '''
    Keywords=['Res','Types']
    List=[]
    def __init__(self, my_file,start_=None,end_=None,type='seq'):
        '''
        Constructor
        '''
        self.Start=0
        self.End=1e24
        
        if start_ != None:
            self.Start=start_
        if end_ != None:
            self.End=end_
        self.my_file=my_file
        self.filename=my_file.name
        self.type=type
        self.N=0
        if self.type =='seq':
            self.iterate=json.loads(self.my_file.readline())
            self.root=self.iterate
            self.beforeStart=True
        else:
            self.root=json.load(self.my_file)
            self.dict=True
            List=[L for L in self.root if L not in JSONIterator.Keywords ]
            self.List=[L for L in List if float(L) >= self.Start and float(L) <= self.End]
            myStart=List[0]
            myEnd=List[len(List)-1]
            if not len(self.List):
                print('Trajectory bounds wrong: Found %-s -- %-s, while given %10.3f -- %10.3f '%(myStart,myEnd,self.Start,self.End))
                sys.exit(1)
            if self.Start < float(myStart):
                self.Start=float(myStart)
            if self.End > float(myEnd):
                self.End=float(myEnd)
            
            self.n=0
            self.end=len(self.List)
            if end_ != None:
                print('   Iterating from step %-11.2f through %-11.2f ' %(self.Start,self.End))
            if start_ == None and end_ == None:
                print('   Iterating over the entire trajectory from step %-11.2f through %-11.2f ' %(self.Start,self.End))
            self.iterate=iter(self.List)

    def __getitem__(self,arg):
        return self.root[arg]

    def __iter__(self):
        return self
    
    def __next__(self):
        if self.type == 'seq':
            line=self.my_file.readline()
            if line:
                self.iterate=json.loads(line)
                self.root=self.iterate
                timeC=next(iter(self.iterate))

                if float(timeC) > self.End:
                    sys.stdout.write('\n')        
                    return None
                
                while self.beforeStart:
                    if float(timeC) > self.Start:
                        self.beforeStart=False
                        print('Start at time = %s ' % timeC)
                        break
                    line=self.my_file.readline()
                    self.iterate=json.loads(line)
                    self.root=self.iterate
                    timeC=next(iter(self.iterate))

                if self.N != 0 and self.N%100 == 0:
                    sys.stdout.write('%10.2f-' % float(timeC))
                    sys.stdout.flush()
                self.N+=1                
                return timeC
            else:
                sys.stdout.write('\n')        
                return None
        else:
            try:
                timeC=next(self.iterate)
                if self.N != 0 and self.N%100 == 0:
                    sys.stdout.write('%10.2f-' % float(timeC))
                    sys.stdout.flush()
                self.N+=1                
                return timeC
            except StopIteration:
                sys.stdout.write('\n')        
                return None
