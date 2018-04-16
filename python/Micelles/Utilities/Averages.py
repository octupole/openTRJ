'''
Created on Dec 13, 2017

@author: marchi
'''
from math import  *
def Wei(x): return 1

class Averages(object):
    '''
    classdocs
    '''

    avg=0.0
    avgW=0.0
    
    def __init__(self, func=Wei):
        '''
        Constructor
        '''
        self.weight=func
        
    def append(self,value,wei=None):
        if not wei:
            w=self.weight(float(value))
            self.avg+=float(value)*w
            self.avgW+=w
        else:
            self.avg+=float(value)*wei
            self.avgW+=wei
            
        
    def __call__(self):
        return self.avg/self.avgW
        
        
if __name__ == "__main__":
    def myW(x):return x*x*x
    A=Averages(myW)
    A.append(3)
    A.append(4.5)
    A.append(6.5)
    A.append(4.5)
    A.append(4.5)
    print(A.average())
    