/*
 * Spline1DInterpolant.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: marchi
 */

#include "Spline1DInterpolant.h"

namespace Spline1D {
Spline1DInterpolant::Spline1DInterpolant(vector<double> xx, vector<double> yy){
	size_t mm{xx.size()};
	x.setlength(mm);
	y.setlength(mm);

	for(int o=0;o<mm;o++){
		x[o]=xx[o];
		y[o]=yy[o];
	}
	spline1dbuildakima(x, y, s);
	Spline=std::bind(spline1dcalc,s,_1);
}
Spline1DInterpolant::Spline1DInterpolant(const vector<vector<double> > & yy){
	size_t mm{yy.size()};
	x.setlength(mm);
	y.setlength(mm);

	int o{0};
	for(auto opset: yy){
		x[o]=opset[0];
		y[o]=opset[1];
		o++;
	}
	spline1dbuildakima(x, y, s);
	Spline=std::bind(spline1dcalc,s,_1);
}
double Spline1DInterpolant::operator()(double x){
	return Spline(x);
}

Spline1DInterpolant::~Spline1DInterpolant() {
	// TODO Auto-generated destructor stub
}

} /* namespace Spline1D */
