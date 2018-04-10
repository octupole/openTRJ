/*
 * SaxsHistograms.cpp
 *
 *  Created on: Feb 19, 2016
 *      Author: marchi
 */

#include "SaxsHistogram.h"



SaxsHistogram::~SaxsHistogram() {
	// TODO Auto-generated destructor stub
}

SaxsHistogram & SaxsHistogram::operator=(SaxsHistogram y){
	hist=y.hist;
	dx=y.dx;
	cut=y.cut;
	Label=y.Label;
	HisX=y.HisX;
	myUnits=y.myUnits;
	myDq=y.myDq;
	return *this;
}
void SaxsHistogram::WriteIt(ostream & fout){
	fout << "# I(Q) Ntot = "<< Label <<endl;
	for(int o=0;o< HisX-1;o++){
		double ddx=dx*static_cast<double>(o);
		if((*this)[o].isZero()) continue;
		double f=(*this)[o].Ratio();
		if(f == 0) continue;
		fout << fixed << setw(8) << setprecision(4) << ddx*myUnits;
		fout << fixed << setw(18) << right << scientific << setprecision(10) << f;
		fout << endl;
	}
}
