/*
 * RhoHistogram.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: marchi
 */

#include "RhoHistogram.h"

namespace Properties {

double RhoHistogram::dx=0.005;
double RhoHistogram::rcut=1.0;
string RhoHistogram::fileout="undefined";

void RhoHistogram::cPrint(ostream & fout){
	Averages();
	for(size_t o{0};o<Histo.size(); o++){
		fout << o*dx*conv_Angstroms<< " " << Avg[o]<<endl;
	}
}
RhoHistogram::RhoHistogram(double R) {
	Rc=R+rcut;
	Histo=vector<double>(static_cast<size_t>(floor(Rc/dx)+1),0.0);
	Avg=Histo;
}
void RhoHistogram::operator ()(double R, double Value){
	if(R > Rc) return;
	if(nHisto.empty()) nHisto=vector<size_t>(static_cast<size_t>(floor(Rc/dx)+1),0);
	size_t n=static_cast<size_t> (floor(R/dx));
	Histo[n]+=Value;
	nHisto[n]++;

}
RhoHistogram & RhoHistogram::operator++(){
	++nstep;
	return *this;
}
RhoHistogram RhoHistogram::operator++(int){
	RhoHistogram temp=*this;
	++*this;
	return temp;
}

void RhoHistogram::operator ()(double R){
	if(R > Rc) return;
	size_t n=static_cast<size_t> (floor(R/dx));
	Histo[n]+=1.0;
}
void RhoHistogram::Averages(double mass){
	this->mass=mass;
	Averages();
}
void RhoHistogram::Averages(){
	auto fact=(1.0_cm*1.0_cm*1.0_cm)/(avogad*1.0_nm*1.0_nm*1.0_nm);
	for(size_t o{0};o<Histo.size();o++){
		auto r0=o*dx;
		auto r1=(o+1)*dx;
		auto vol0=4.0*M_PI*r0*r0*r0/3.0;
		auto vol1=4.0*M_PI*r1*r1*r1/3.0;
		auto dV=vol1-vol0;
		double rho=(Histo[o]*mass)/dV;
		rho*=fact.val;
		Avg[o]=rho/(double) nstep;
	}
}
vector<double> RhoHistogram::operator[](size_t n){
	vector<double> out{0.0,0.0};

	out[0]=n*dx;
	out[1]=Avg[n];
	return out;
}

RhoHistogram::~RhoHistogram() {
	// TODO Auto-generated destructor stub
}

} /* namespace fftwpp */
