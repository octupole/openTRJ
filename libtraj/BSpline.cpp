/*
 * BSpline.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: marchi
 */

#include "BSpline.h"
int BSpline::order=4;

void BSpline::allocate(){
	theta.x=vector<double>(order);
	theta.dx=vector<double>(order);
}
spline & BSpline::operator ()(const double w){
	Init(w);
	for(auto o=2;o<order-1;o++)
		OnePass(w,o);
	Diff();
	OnePass(w,order-1);
	return theta;
}

void BSpline::Init(const double w){
	for(auto i=2;i<order;i++) theta.x[i]=0.0;
	theta.x[1]=w;
	theta.x[0]=1.0-w;
};
void BSpline::OnePass(const double w ,const int k1){
	int k=k1+1;
	double div=1.0/static_cast<double>(k-1);
	theta.x[k-1]=div*w*theta.x[k-2];
	for(auto j=1;j<=k-2;j++)
		theta.x[k-j-1]=div*((w+j)*theta.x[k-j-2]+(k-j-w)*theta.x[k-j-1]);
	theta.x[0]=div*(1.0-w)*theta.x[0];

/*
	double div;
	int K1=k+1;
	div=1.0/ static_cast<double>(K1-1);
	theta.x[k]=div*w*theta.x[k-1];
	for(auto o=0;o<k-2;o++){
		theta.x[k-o-1]=div*((w+o+1)*theta.x[k-o-2]+(k-o-w)*theta.x[k-o-1]);
	}
	theta.x[0]*=div*(1.0-w)*theta.x[0];
	*/
}
void BSpline::Diff(){
	theta.dx[0]=-theta.x[0];
	for(auto o=1;o<order;o++)
		theta.dx[o]=theta.x[o-1]-theta.x[o];
}


BSpline::BSpline() {
	// TODO Auto-generated constructor stub
	allocate();
}

BSpline::BSpline(int myorder) {
	// TODO Auto-generated constructor stub
	order=myorder;
	allocate();
}

BSpline::~BSpline() {
	// TODO Auto-generated destructor stub
}

