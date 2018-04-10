/*
 * IteratorMAtoms.cpp
 *
 *  Created on: Dec 17, 2015
 *      Author: marchi
 */

#include "IteratorVoronoi.h"

namespace myiterators {

IteratorVoronoi::IteratorVoronoi(pointer x, ifstream * y,long int nstrt,long int nnd): p{x},
finx{y},nstart{nstrt}, nend{nnd}{
	ios::streampos whereIwas=finx->tellg();
	finx->seekg(0,finx->end);
	len=finx->tellg();
	finx->seekg(0,finx->beg);
	finx->seekg(whereIwas,finx->beg);
}


IteratorVoronoi & IteratorVoronoi::operator++(){
	ifstream & fin=*finx;
	while(finx->tellg() < len && !finx->fail()){
		if(ntime < nstart) {fin>> *p;ntime++;;continue;}
		if(nend != -1 && ntime > nend) break;

		fin >> *p;
		ntime++;
		return *this;
	}
	(*this)=IteratorVoronoi();
	return *this;
}



} /* namespace myiterators */
