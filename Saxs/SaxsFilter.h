/*
 * SaxsFilter.h
 *
 *  Created on: Feb 24, 2017
 *      Author: marchi
 */

#ifndef SRC_SAXSFILTER_H_
#define SRC_SAXSFILTER_H_
#include "Array.h"
#include <iostream>
#include <string>
#include <complex>

using std::cout;
using std::endl;
using std::string;
using Complex=std::complex<double>;

using namespace Array;

class SaxsFilter {
	array3<Complex> * Filter_k{nullptr};
	size_t nnx{0},nny{0},nnz{0};
	void __Allocate();
public:
	SaxsFilter();
	SaxsFilter(size_t,size_t,size_t);
	array3<Complex> & operator*(array3<Complex> & );
	void Allocate(size_t,size_t,size_t);
	bool allocated(){return Filter_k;}
	virtual ~SaxsFilter();
};

#endif /* SRC_SAXSFILTER_H_ */
