/*
 * AtomIndex.h
 *
 *  Created on: May 26, 2011
 *      Author: marchi
 */

#ifndef ATOMINDEX_H_
#define ATOMINDEX_H_
#include <cstdlib>
#include <iostream>
#include "Ftypedefs.h"

using namespace Typedefs;
class AtomIndex {
	int * idx;
	int nidx;
public:
	AtomIndex():idx(NULL),nidx(0){};
	AtomIndex(int nd):nidx(nd){idx=new int [nidx];};
	AtomIndex(AtomIndex & cind){
		nidx=cind.nidx;
		idx=new int[nidx];
		for(int i=0;i<nidx;i++) idx[i]=cind[i];
	}
	AtomIndex(int nd, int * cind):nidx(nd){
		idx=new int [nidx];
		for(int i=0;i<nidx;i++) idx[i]=cind[i];
	};
	AtomIndex & operator=(AtomIndex & c){
		if(nidx) delete [] idx;
		nidx=c.nidx;
		idx=new int[nidx];
		for(int i=0;i<nidx;i++) idx[i]=c[i];
		return *this;
	}
	const int * getIdx() const {return idx;}
	virtual ~AtomIndex(){if(idx) delete [] idx;};
	void set(int nd, const int * cind){
		nidx=nd;
		if(idx) delete [] idx;
		idx=new int [nidx];
		for(int i=0;i<nidx;i++) idx[i]=cind[i];
	};
	void Mol(const int n,const t_blockx & mol){
		nidx=mol.index[n+1]-mol.index[n];
		if(idx) delete [] idx;
		idx=new int [nidx];
		int m=0;
		for(int i=mol.index[n];i<mol.index[n+1];i++){
			idx[m++]=i;
		}
	}
	void MolCompl(const int n,const t_blockx & mol){
		nidx=mol.index[mol.nr]-mol.index[0]-(mol.index[n+1]-mol.index[n]);
		if(idx) delete [] idx;
		idx=new int [nidx];
		int j=0;
		for(int m=0;m<mol.nr;m++){
			if(m==n) continue;
			for(int i=mol.index[m];i<mol.index[m+1];i++)
				idx[j++]=i;
		}
	}
	void MolTot(const t_blockx & mol){
		nidx=mol.index[mol.nr]-mol.index[0];
		if(idx) delete [] idx;
		idx=new int [nidx];
		for(int i=mol.index[0];i<mol.index[mol.nr];i++)
			idx[i]=i;
	}
	void Print(){
		for(int i=0;i<nidx;i++)
			std::cout << idx[i]<< std::endl;
	}
	int getN() const {return nidx;};
	void setN(const int n){nidx=n;};
	int & operator[](const int i) const {return idx[i];};
	int & getI(const int i) const {return idx[i];};

};

#endif /* ATOMINDEX_H_ */
