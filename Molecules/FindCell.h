/*
 * FindCell.h
 *
 *  Created on: May 18, 2016
 *      Author: marchi
 */

#ifndef TRJLIB_FINDCELL_H_
#define TRJLIB_FINDCELL_H_
#include "MyUtilClass.h"
#include <vector>
#include <random>
#include <cmath>
#include "Ftypedefs.h"
#include <algorithm>
#include <iostream>

using namespace MATRIX;
using namespace DVECT;
using std::vector;
using std::cout;
using std::endl;

template <typename T>
class FindCell {
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	const int Maximum{10000};
	MMatrix<T> CO{0.0};
public:
	FindCell(MMatrix<T> co):CO{co}{};
	Dvect Run(const vector<vector<int> > &, const vector<vector<int> > &,const vector<Dvect> , int);
	virtual ~FindCell();
};

#endif /* TRJLIB_FINDCELL_H_ */
