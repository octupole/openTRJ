/*
 * SaxsHistogramSpline.h
 *
 *  Created on: Mar 21, 2016
 *      Author: marchi
 */

#ifndef TRJLIB_SAXSHISTOGRAMSPLINE_H_
#define TRJLIB_SAXSHISTOGRAMSPLINE_H_

#include "SaxsHistogram.h"
#include "Spline1DInterpolant.h"
/** \brief Derived from SaxsHistogram uses splines instead of linear interpolation
 *
 */
class SaxsHistogramSpline: public SaxsHistogram {
	virtual void WriteIt(ostream &);
public:
	SaxsHistogramSpline():SaxsHistogram(){};
	SaxsHistogramSpline(double a,double b, double c):SaxsHistogram(a,b,c){};
	virtual ~SaxsHistogramSpline();
};

#endif /* TRJLIB_SAXSHISTOGRAMSPLINE_H_ */
