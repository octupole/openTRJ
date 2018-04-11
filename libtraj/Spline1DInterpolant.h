/*
 * Spline1DInterpolant.h
 *
 *  Created on: Mar 22, 2016
 *      Author: marchi
 */

#ifndef TRJLIB_SPLINE1DINTERPOLANT_H_
#define TRJLIB_SPLINE1DINTERPOLANT_H_
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <set>
#include "stdafx.h"
#include "interpolation.h"
using namespace alglib;
using namespace std::placeholders;
using std::bind;

using std::function;
using std::vector;
using std::set;
namespace Spline1D {

class Spline1DInterpolant {
	real_1d_array x,y;
	spline1dinterpolant s;
public:
	Spline1DInterpolant()=delete;
	Spline1DInterpolant(vector<double>,vector<double>);
	Spline1DInterpolant(const vector<vector<double> >& );
	double operator()(double);
	function<double(const double)> Spline;
	virtual ~Spline1DInterpolant();
};

} /* namespace Spline1D */

#endif /* TRJLIB_SPLINE1DINTERPOLANT_H_ */
