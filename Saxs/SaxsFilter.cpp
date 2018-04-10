/*
 * SaxsFilter.cpp
 *
 *  Created on: Feb 24, 2017
 *      Author: marchi
 */

#include "SaxsFilter.h"
void SaxsFilter::__Allocate(){
	Filter_k=new array3<Complex>(nnx,nny,nnz);
}

SaxsFilter::SaxsFilter(size_t mx,size_t my,size_t mz):nnx{mx},nny{my},nnz{mz}{
	__Allocate();
}
SaxsFilter::SaxsFilter() {}

void SaxsFilter::Allocate(size_t mx,size_t my,size_t mz){
	nnx=mx;nny=my;nnz=mz;
	__Allocate();
}
array3<Complex> & SaxsFilter::operator*(array3<Complex> & x){
	if(!Filter_k) return x;
	try{
		if(x.Size() != Filter_k->Size()) throw "Trying to apply a Saxs filter but dimensions are wrong.";
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
	for(auto o=0;o<nnx;o++)
		for(auto p=0;p<nny;p++)
			for(auto q=0;q<nnz;q++)
				x[o][p][q]*=(*Filter_k)[o][p][q];
	return x;
}

SaxsFilter::~SaxsFilter() {
	if(Filter_k) {
		Filter_k->Deallocate();
		delete Filter_k;
	}
}

