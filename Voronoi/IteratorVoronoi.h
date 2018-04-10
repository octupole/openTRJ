/*
 * IteratorMAtoms.h
 *
 *  Created on: Dec 17, 2015
 *      Author: marchi
 */

#ifndef SRC_ITERATORVORONOI_H_
#define SRC_ITERATORVORONOI_H_

#include <iterator>
#include <fstream>
#include "Voronoi.h"
#include "NewMPI.h"
using namespace Voro;
using std::ifstream;
namespace myiterators {
/** \brief Iterator class to iterate on the trajectory steps
 *
 */

class IteratorVoronoi: public std::iterator<std::output_iterator_tag, Voronoi >{
protected:
	using pointer=Voronoi *;
	using reference=Voronoi &;
	pointer p{nullptr};
	ifstream * finx{nullptr};
	long int ntime{0};
	long int nstart{0},nend{-1};
	vector<long int> Nst,Nnd;
	ios::streampos len;
public:
	IteratorVoronoi(){};
	IteratorVoronoi(pointer , ifstream *,long int ,long int );
	IteratorVoronoi(const IteratorVoronoi & mit): p{mit.p},finx{mit.finx}, nstart{mit.nstart}
		, nend{mit.nend}, len{mit.len}{};
	virtual void Initialize(Parallel::NewMPI * y){};
	virtual IteratorVoronoi & operator++();
	bool isReferenced(){return p;}
	long int getTime(){return ntime-1;}

	virtual IteratorVoronoi operator++(int) {IteratorVoronoi tmp(*this); operator++(); return tmp;}

	reference operator*() {return *p;}
	pointer  operator()() {return p;}
	virtual ~IteratorVoronoi(){
		delete p;
	}
};
} /* namespace myiterators */

#endif /* SRC_ITERATORVORONOI_H_ */
