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
#include "stdafx.h"
#include "interpolation.h"
#include "histograms.hpp"
using std::bind;
using std::function;
namespace Spline1D {
/** \brief Interpolates the I(q)
 *
 */
class Spline1DInterpolant {
	alglib::real_1d_array x,y;
	alglib::spline1dinterpolant s;
	double Dq{0},myUnits{1.0},cutoff{0};
public:
	Spline1DInterpolant()=delete;
	Spline1DInterpolant(Histogram1D *, double,double);
	double lowLimit()const {return x[0];}
	double operator()(double);
	Spline1DInterpolant & operator-=(const Spline1DInterpolant &);
	function<double(const double)> Spline;
	virtual ~Spline1DInterpolant();
	friend ostream & operator<<(ofstream &,Spline1DInterpolant &);
};

} /* namespace Spline1D */

#endif /* TRJLIB_SPLINE1DINTERPOLANT_H_ */
