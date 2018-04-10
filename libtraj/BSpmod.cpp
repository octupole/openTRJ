/*
 * BSpmod.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: marchi
 */

#include "BSpmod.h"

size_t BSpmod::nx=0;
size_t BSpmod::ny=0;
size_t BSpmod::nz=0;
size_t BSpmod::ndim=0;
size_t BSpmod::order=0;

const double twopi=M_PI*2.0;
const double tiny=1.0e-7;
void BSpmod::Inverse(){
	for(auto & v: BSp.x){
		v=1.0/v;
	}
	for(auto & v: BSp.y){
		v=1.0/v;
	}
	for(auto & v: BSp.z){
		v=1.0/v;
	}

}
void BSpmod::load_moduli(){
	BSpline tmp;
	double w=0.0;
	spline d=BSpline{}(w);
	vector<double> A0(ndim,0.0);
	for(auto o=1;o<order+1;o++)
		A0[o]=d.x[o-1];
	BSp.x=DFTmod(A0,nx);
	Gamma(BSp.x);
	BSp.y=DFTmod(A0,ny);
	Gamma(BSp.y);
	BSp.z=DFTmod(A0,nz);
	Gamma(BSp.z);
	Inverse();
}
vector<double> BSpmod::DFTmod(const vector<double> & A, size_t Ndim){
	vector<double> bsp(Ndim,0.0);
	for(auto o=0;o<Ndim;o++){
		double sum1=0.0,sum2=0.0;
		for(auto p=0;p<order+1;p++){
			double arg=twopi*static_cast<double>(o*p)/static_cast<double>(Ndim);
			sum1+=A[p]*cos(arg);
			sum2+=A[p]*sin(arg);
		}
		bsp[o]=sum1*sum1+sum2*sum2;
	}
	for(auto o=1;o<Ndim-1;o++)
		bsp[o]=bsp[o]>tiny?bsp[o]:0.5*(bsp[o-1]+bsp[o+1]);
	return bsp;
}

void BSpmod::Gamma(vector<double> & bsp){

	size_t Ndim=bsp.size();
	size_t nf=Ndim/2;

	auto GammaSum=[&Ndim](const int m, const int order)->double{
		if(!m) return 1.0;
		double gsum=1.0;
		double x=M_PI*static_cast<double>(m)/static_cast<double>(Ndim);
		for(auto k=1;k<=KCUT;k++){
			gsum+=pow(x/(x+M_PI*k),order)+pow(x/(x-M_PI*k),order);
		}
		return gsum;
	};

	int order2=2*order;
	for(auto k=0;k<Ndim;k++){
		double lambda=1.0;
		int m=k<nf?k:k-Ndim;
		if(m != 0) {
			double gsum=GammaSum(m,order);
			double gsum2=GammaSum(m,order2);
			lambda=gsum/gsum2;
		}

		bsp[k]/=(lambda*lambda);
	}
}

BSpmod::~BSpmod() {
	// TODO Auto-generated destructor stub
}

