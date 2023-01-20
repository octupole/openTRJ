/*
 * IteratorMAtoms.h
 *
 *  Created on: Dec 17, 2015
 *      Author: marchi
 */

#ifndef SRC_ITERATORMATOMS_H_
#define SRC_ITERATORMATOMS_H_

#include <iterator>
#include "MAtoms.h"
#include "FstreamC.h"
#include "FstreamF.h"
#include "NewMPI.h"

namespace myiterators {
/** \brief Iterator class to iterate on the trajectory steps
 *
 */
class IteratorMAtoms: public std::iterator<std::output_iterator_tag, MAtoms>{
protected:
	pointer p{nullptr};
	Fstream * finx{nullptr};
	long int ntime{0};
	long int nstart{0},nend{-1},nskip{1};
	vector<long int> Nst,Nnd,Nskp;
	streampos len;
public:
	IteratorMAtoms(){};
	IteratorMAtoms(pointer x, Fstream * y,long int nstrt,long int nnd,long int nskp): p{x}
	,finx{y}, nstart{nstrt}, nend{nnd}, nskip{nskp}{
		finx->seekg(0,"end");
		len=finx->tellg();
		finx->Rewind();
	};
	IteratorMAtoms(const IteratorMAtoms & mit): p{mit.p},finx{mit.finx}, nstart{mit.nstart}
		, nend{mit.nend}, nskip{mit.nskip}, len{mit.len}{};
	virtual void Initialize(Parallel::NewMPI * y){};
	virtual IteratorMAtoms & operator++(){
		Fstream & fin=*finx;
		while(finx->tellg() < len){
			if(ntime < nstart) {fin+=*p;ntime++;;continue;}
			if(nend != -1 && ntime > nend) break;
			if(ntime == 0) {
				fin+=*p;
				ntime++;
				continue;
			}

			if(!((ntime-nstart)%nskip)) {
				fin >> *p;

				ntime++;
				return *this;
			} else {
				fin+=*p;
			}
			ntime++;
		}
		(*this)=IteratorMAtoms();
		return *this;
	}
	bool isReferenced(){return p;}
	long int getTime(){return ntime-1;}

	virtual IteratorMAtoms operator++(int) {IteratorMAtoms tmp(*this); operator++(); return tmp;}

	reference operator*() {return *p;}
	pointer  operator()() {return p;}
	virtual ~IteratorMAtoms(){
		delete p;
	}
};
} /* namespace myiterators */

#endif /* SRC_ITERATORMATOMS_H_ */
