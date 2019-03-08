/*
 * RhoHistogram.cpp
 *
 *  Created on: Mar 7, 2019
 *      Author: marchi
 */

#include "RhoHistogram.h"

namespace Properties {

string RhoHistogram::headerXVG=R"(
# Grace project file
#
@version 50122
@page size 792, 612
@page scroll 5%
@page inout 5%
@link page off
@map font 13 to "ZapfDingbats", "ZapfDingbats"
@map font 4 to "Helvetica", "Helvetica"
@map font 6 to "Helvetica-Bold", "Helvetica-Bold"
@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"
@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"
@map font 0 to "Times-Roman", "Times-Roman"
@map font 2 to "Times-Bold", "Times-Bold"
@map font 1 to "Times-Italic", "Times-Italic"
@map font 3 to "Times-BoldItalic", "Times-BoldItalic"
@map font 8 to "Courier", "Courier"
@map font 10 to "Courier-Bold", "Courier-Bold"
@map font 9 to "Courier-Oblique", "Courier-Oblique"
@map font 12 to "Symbol", "Symbol"
@map color 0 to (255, 255, 255), "white"
@map color 1 to (0, 0, 0), "black"
@map color 2 to (255, 0, 0), "red"
@map color 3 to (0, 255, 0), "green"
@map color 4 to (0, 0, 255), "blue"
@map color 5 to (255, 255, 0), "yellow"
@map color 6 to (188, 143, 143), "brown"
@map color 7 to (220, 220, 220), "grey"
@map color 8 to (148, 0, 211), "violet"
@map color 9 to (0, 255, 255), "cyan"
@map color 10 to (255, 0, 255), "magenta"
@map color 11 to (255, 165, 0), "orange"
@map color 12 to (114, 33, 188), "indigo"
@map color 13 to (103, 7, 72), "maroon"
@map color 14 to (64, 224, 208), "turquoise"
@map color 15 to (0, 139, 0), "green4"
@reference date 0
@date wrap off
@date wrap year 1950
@default linewidth 1.0
@default linestyle 1
@default color 1
@default pattern 1
@default font 0
@default char size 1.000000
@default symbol size 1.000000
@default sformat "%.8g"
@background color 0
@page background fill on
@timestamp off
@timestamp 0.03, 0.03
@timestamp color 1
@timestamp rot 0
@timestamp font 0
@timestamp char size 1.000000
@timestamp def "Mon Oct  8 15:33:20 2012"
)";
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
void RhoHistogram::Reduce(Parallel::NewMPI * y){
	if(!y->Get_Size()) return;
	int MPIsize=y->Get_Size();
	int * i_recv=new int [MPIsize];

	int nstep=this->nstep;
	y->ReduceSum(&nstep,1);
	this->nstep=nstep;

	int dim=nHisto.size();
	y->AllGather(1,&dim,i_recv);
	dim=*std::max_element(i_recv,i_recv+MPIsize);
	nHisto.resize(dim);
	auto ip0=&nHisto[0];
	y->ReduceSum(ip0,dim);

	dim=Histo.size();
	y->AllGather(1,&dim,i_recv);
	dim=*std::max_element(i_recv,i_recv+MPIsize);
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
	fout << RhoHistogram::headerXVG<<endl;
	for(auto it{y.begin()};it != y.end();it++){
		string key=it->first;
		fout << *(y[key]) ;
		fout << "&"<<endl;
	}
	return fout;
};

} /* namespace fftwpp */
