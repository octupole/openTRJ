/*
 * SaxsHistogramSpline.cpp
 *
 *  Created on: Mar 21, 2016
 *      Author: marchi
 */

#include "SaxsHistogramSpline.h"


void SaxsHistogramSpline::WriteIt(ostream & fout){
	Spline1D::Spline1DInterpolant spline0(this,dx,myUnits);
	fout << "# I(Q) Ntot = "<< Label <<endl;
	for(int o=0;o< (HisX-1)*dx/myDq+1;o++){
		double ddx=myDq*static_cast<double>(o);
		double f{spline0(ddx)};
		if(f == 0) continue;
		fout << fixed << setw(8) << setprecision(4) << ddx*myUnits;
		fout << fixed << setw(18) << right << scientific << setprecision(10) << f;
		fout << endl;
	}
}


SaxsHistogramSpline::~SaxsHistogramSpline() {
	// TODO Auto-generated destructor stub
}

