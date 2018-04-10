/*
 * Spline1DInterpolant.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: marchi
 */

#include "Spline1DInterpolant.h"
using namespace alglib;
using namespace std::placeholders;

namespace Spline1D {
function<double(real_1d_array,real_1d_array)> LeastSquare=[](real_1d_array x0, real_1d_array y0)->double{
	const ae_int_t Dims=7;
	real_1d_array xx;
	real_1d_array yy;
	xx.setlength(Dims);
	yy.setlength(Dims);

	ae_int_t info,m{2};
	barycentricinterpolant p;
	polynomialfitreport rep;

	int op{0},pp{1};
	vector<double> X,Y;

	while(op <= Dims){
	  double x{x0[pp]},y{y0[pp]};
	  pp++;
	  if(y != 0.0) {
	    X.push_back(x*x);
	    Y.push_back(log(y));
	    op++;
	  }
	}
	for(auto o=0;o<Dims;o++){
	  xx[o]=X[o];
	  yy[o]=Y[o];
	}

	polynomialfit(xx, yy, Dims,m, info, p, rep);
	double t=barycentriccalc(p,0);

	return exp(t);
};
Spline1DInterpolant::Spline1DInterpolant(Histogram1D * qdfx,double dq, double units):Dq{dq}, myUnits{units} {
	cutoff=qdfx->dx*(qdfx->HisX-1);
	int mm{0};
	for(int o=1;o< qdfx->HisX-1;o++){
		if((*qdfx)[o].isZero()) continue;
		if((*qdfx)[o].Ratio() == 0) continue;
		mm++;
	}
	x.setlength(mm);
	y.setlength(mm);
	mm=0;
	for(int o=1;o< qdfx->HisX-1;o++){
		if((*qdfx)[o].isZero()) continue;
		if((*qdfx)[o].Ratio() == 0) continue;

		x[mm]=qdfx->dx*static_cast<double>(o);
		y[mm]=(*qdfx)[o].Ratio();
		mm++;
	}
	spline1dbuildakima(x, y, s);
	Spline=std::bind(spline1dcalc,s,_1);
}
double Spline1DInterpolant::operator()(double x){
	return Spline(x);
}
Spline1DInterpolant & Spline1DInterpolant::operator-=(const Spline1DInterpolant & y0){
	double lowL=std::max(x[0],y0.lowLimit());
	int mm{0};
	real_1d_array w,z;
	int length=(int)cutoff/Dq+1;
	w.setlength(length);
	z.setlength(length);
	z[0]=0;
	w[0]=0;
	for(auto o=1;o<length;o++){
		z[o]=o*Dq;
		w[o]=Spline(z[o])-y0.Spline(z[o]);
	}
	w[0]=LeastSquare(z,w);
	x=z;y=w;
	spline1dbuildakima(x, y, s);
	Spline=std::bind(spline1dcalc,s,_1);
	return *this;
}

ostream & operator<<(ofstream & fout,Spline1DInterpolant & y){
	fout << "# Interpolated SAXS data " << endl;
	int length=y.cutoff/y.Dq;
	for(int o=0;o< length;o++){
		double x=y.Dq*static_cast<double>(o);
		double f=y.Spline(x);
		if(!f) continue;
		fout << fixed << setw(8) << setprecision(5) << x*y.myUnits;
		fout << fixed << setw(12) << right << scientific << setprecision(4) << f;
		fout << endl;
	}
	return fout;
}
Spline1DInterpolant::~Spline1DInterpolant() {
	// TODO Auto-generated destructor stub
}

} /* namespace Spline1D */
