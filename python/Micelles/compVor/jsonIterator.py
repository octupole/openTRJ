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
            self.lines=iter(self.my_file.readlines())
            self.iterate=json.load(next(self.lines))
            self.root=self.iterate

    def __getitem__(self,arg):
        return self.root[arg]

    def __iter__(self):
        return self
    
    def __next__(self):
        if self.type == 'seq':
            line=self.my_file.readline()
        else:
            line=next(self.lines)

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