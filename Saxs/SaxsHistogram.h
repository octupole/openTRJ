/*
 * SaxsHistograms.h
 *
 *  Created on: Feb 19, 2016
 *      Author: marchi
 */

#ifndef SRC_SAXSHISTOGRAM_H_
#define SRC_SAXSHISTOGRAM_H_

#include "histograms.hpp"
#include "stdafx.h"
#include "interpolation.h"

using namespace alglib;
/** Derived class of Histogram1D specific for saxs, where units are inverse of distances
 *
 */
class SaxsHistogram: public Histogram1D{
protected:
	double myUnits{0.0};
	double myDq{0.0};
	virtual void WriteIt(ostream &);
public:
	SaxsHistogram():Histogram1D(){};
	SaxsHistogram(double a,double b, double units):myUnits{units},Histogram1D(a,b){};
	virtual SaxsHistogram & operator^=(const SaxsHistogram &);
	SaxsHistogram & operator=(SaxsHistogram );
	double gmyUnits(){return myUnits;}
	void setDq(double dq){myDq=dq;}
	virtual ~SaxsHistogram();
	friend ostream & operator<<(ostream & fout, SaxsHistogram & y){
		y.WriteIt(fout);
		return fout;
	}

};

#endif /* SRC_SAXSHISTOGRAM_H_ */
