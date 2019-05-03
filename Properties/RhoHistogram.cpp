/*
 * RhoHistogram.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: marchi
 */

#include "RhoHistogram.h"

namespace Properties {


double RhoHistogram::dx=0.005;
double RhoHistogram::rcut=-1.0;
string RhoHistogram::fileout="undefined";

void RhoHistogram::cPrint(ostream & fout){
	Averages();
	for(size_t o{0};o<Histo.size(); o++){
		fout << o*dx*conv_Angstroms<< " " << Avg[o]<<endl;
	}
}
RhoHistogram::RhoHistogram(double R) {
	Rc=rcut<0?R:rcut;
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
void RhoHistogram::Reduce(Parallel::NewMPI * y){
	if(!y->Get_Size()) return;
	int MPIsize=y->Get_Size();
	int * i_recv=new int [MPIsize];

	int nstep=this->nstep;
	y->ReduceSum(&nstep,1);
	this->nstep=nstep;

	int dim=nHisto.size();
	y->AllGather(1,&dim,i_recv);
	dim=*std::min_element(i_recv,i_recv+MPIsize);
	nHisto.resize(dim);
	auto ip0=&nHisto[0];
	y->ReduceSum(ip0,dim);

	dim=Histo.size();
	y->AllGather(1,&dim,i_recv);
	dim=*std::min_element(i_recv,i_recv+MPIsize);
	Histo.resize(dim);
	Avg.resize(dim);
	auto ip1=&Histo[0];
	y->ReduceSum(ip1,dim);
}
RhoHistogram::~RhoHistogram() {
	// TODO Auto-generated destructor stub
}
ostream & operator << (ostream & fout , RhoHistogram  & y ){
	y.cPrint(fout);
	return fout;
}
ostream & operator << (ostream & fout , map<string,RhoHistogram *> & y ){
	int N{0};
	std::stringstream ss;
	for(auto it=y.begin();it != y.end();it++){
		string key=it->first;
		ss << "@    s"<<N<<" legend  \""<<key<<"\""<<endl;
		N++;
	}
	fout << ss.str();
	for(auto it=y.begin();it != y.end();it++){
		string key=it->first;
		fout << *(y[key]) ;
		fout << "&"<<endl;
	}
	return fout;
};

} /* namespace fftwpp */
