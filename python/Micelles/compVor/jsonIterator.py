'''
Created on Dec 17, 2017

@author: marchi
'''
import sys 

class JSONIterator(object):
    '''
    classdocs
    '''
    Keywords=['Res','Types']
    List=[]
    def __init__(self, root,start_=None,end_=None):
        '''
        Constructor
        '''
        self.Start=0
        self.End=1e24
        
        if start_ != None:
            self.Start=start_
        if end_ != None:
            self.End=end_
        self.root=root
        List=[L for L in root if L not in JSONIterator.Keywords ]
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
            
        self.iterate=iter(self.List)
        self.n=0
        self.end=len(self.List)
        if end_ != None:
            print('   Iterating from step %-11.2f through %-11.2f ' %(self.Start,self.End))
        if start_ == None and end_ == None:
            print('   Iterating over the entire trajectory from step %-11.2f through %-11.2f ' %(self.Start,self.End))
            
    def __iter__(self):
        return self
    
    def next(self):
        try:
            return next(self.iterate)
        except StopIteration:
            return None
#             
#         print(self.n,self.end)
#         if self.n < self.end:
#             if self.List[self.n] in self.Keywords:
#                 self.n+=1
#                 self.next()
#             else:
#                 while True:
#                     self.item=self.List[self.n]
#                     self.n+=1
#                     if float(self.item) >= self.Start and float(self.item) <= self.End:
#                         break
#         return self.item
#         