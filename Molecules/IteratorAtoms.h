/*
 * IteratorMAtoms.h
 *
 *  Created on: Dec 17, 2015
 *      Author: marchi
 */

#ifndef SRC_ITERATORMATOMS_H_
#define SRC_ITERATORMATOMS_H_

#include <iterator>
#include "Atoms.h"
#include "FstreamC.h"
#include "FstreamF.h"
#include "NewMPI.h"

namespace myiterators {
/** \brief Iterator class to iterate on the trajectory steps
 *
 */
template <typename T>
class IteratorAtoms: public std::iterator<std::output_iterator_tag, Atoms<T> >{
protected:
	using pointer=Atoms<T> *;
	using reference=Atoms<T> &;
	pointer p{nullptr};
	Fstream * finx{nullptr};
	long int ntime{0};
	long int nstart{0},nend{-1},nskip{1};
	vector<long int> Nst,Nnd,Nskp;
	ios::streampos len;
public:
	IteratorAtoms(){};
	IteratorAtoms(pointer , Fstream *,long int ,long int ,long int);
	IteratorAtoms(const IteratorAtoms & mit): p{mit.p},finx{mit.finx}, nstart{mit.nstart}
		, nend{mit.nend}, nskip{mit.nskip}, len{mit.len}{};
	virtual void Initialize(Parallel::NewMPI * y){};
	virtual IteratorAtoms & operator++();
	bool isReferenced(){return p;}
	long int getTime(){return ntime-1;}

	virtual IteratorAtoms operator++(int) {IteratorAtoms tmp(*this); operator++(); return tmp;}

	reference operator*() {return *p;}
	pointer  operator()() {return p;}
	virtual ~IteratorAtoms(){
		delete p;
	}
};
} /* namespace myiterators */

#endif /* SRC_ITERATORMATOMS_H_ */
