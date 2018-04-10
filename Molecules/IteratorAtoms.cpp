/*
 * IteratorMAtoms.cpp
 *
 *  Created on: Dec 17, 2015
 *      Author: marchi
 */

#include <IteratorAtoms.h>

namespace myiterators {

template <typename T>
IteratorAtoms<T>::IteratorAtoms(pointer x, Fstream * y,long int nstrt,long int nnd,long int nskp): p{x}
,finx{y}, nstart{nstrt}, nend{nnd}, nskip{nskp}{
	finx->seekg(0,"end");
	len=finx->tellg();
	finx->Rewind();
}

template <typename T>
IteratorAtoms<T> & IteratorAtoms<T>::operator++(){
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
	(*this)=IteratorAtoms();
	return *this;
}


template class IteratorAtoms<float>;
template class IteratorAtoms<double>;

} /* namespace myiterators */
