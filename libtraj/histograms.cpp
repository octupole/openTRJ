/*
 * histograms.cpp
 *
 *  Created on: Jun 25, 2015
 *      Author: marchi
 */

#include "../libtraj/histograms.hpp"



Histogram1D::Histogram1D(double dx0,double cut0): dx(dx0), cut(cut0), HisX(0), Label(" "){
	int nh=static_cast<int>(cut/dx)+2;
	HisX=nh;
	hist=vector<hist1D>(2*nh+1);
}
void Histogram1D::operator()(double dx0,double cut0){
	dx=dx0;cut=cut0;
	int nh=static_cast<int>(cut/dx)+2;
	HisX=nh;
	hist=vector<hist1D>(2*nh+1);
}
hist1D & Histogram1D::operator[](size_t n){
	try{
		if(n > hist.size()) throw "Something wrong with histogram ";
		} catch(const char * s){
			cout << s << " " << n << " " << hist.size() <<endl;
			exit(1);
		}
	return hist.at(HisX+n);
}

Histogram1D & Histogram1D::operator*=(double y){
	for(auto & op: hist)
		op*=y;
	return *this;
}
Histogram1D & Histogram1D::operator/=(double y){
	for(auto & op: hist)
		op/=y;
	return *this;
}
Histogram1D & Histogram1D::operator*=(size_t y){
	double y1=static_cast<double>(y);
	for(auto & op: hist)
		op*=y1;
	return *this;
}

Histogram1D & Histogram1D::operator=(const Histogram1D & z){
	hist=z.hist;
	dx=z.dx;
	cut=z.cut;
	Label=z.Label;
	HisX=z.HisX;
	return *this;
}
bool Histogram1D::operator==(const Histogram1D & z){
	return (hist.size() == z.hist.size());
}
Histogram1D & Histogram1D::operator+=(const Histogram1D & z){
	dx=z.dx;
	cut=z.cut;
	Label=z.Label;
	HisX=z.HisX;
	try{
		if(hist.size() == 0){
			hist=vector<hist1D>(z.hist.size(),hist1D());
			} else if(hist.size() != z.hist.size()) throw string("The hist1D dimension are not identical ");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}

	for(size_t o=0; o < hist.size(); o++) hist[o]+=z.hist[o];
	return *this;
}
void Histogram1D::clear(){
	hist.clear();
	dx=0;
	cut=0;
	Label="";
	HisX=0;
}

void Histogram1D::WriteIt(ostream & fout){
	fout << "# Radial density for Selection labeled " + Label + " "<< endl;
	fout << "#  R       Rho(r) " << endl;
	for(size_t o=0;o< HisX-2;o++){
		double ddx=dx*static_cast<double>(o);
		double f=(*this)[o].Ratio();
		if(!f) continue;
		fout << fixed << setw(8) << setprecision(5) << ddx;
		fout << fixed << setw(12) << right << scientific << setprecision(4) << f;
		fout << endl;
	}

}

ostream & operator<<(ostream & fout, Histogram1D & y){
	y.WriteIt(fout);
	return fout;
}

void PairCorr1D::WriteIt(ostream & fout){
	fout << "# Radial density for Selection labeled " + Label + " "<< endl;
	fout << "#  R       Rho(r) " << endl;
	for(size_t h=0;h< HisX-2;h++){
		double r=dx*static_cast<double>(h);
		double vol=4.0*M_PI*r*r*dx;
		double fact=1.0/static_cast<double>(nMol)/vol;
		double f=(*this)[h].Ratio()*fact;
		if(!f) continue;
		fout << fixed << setw(8) << setprecision(5) << r;
		fout << fixed << setw(12) << right << scientific << setprecision(4) << f;
		fout << endl;
	}

}

ostream & operator<<(ostream & fout, PairCorr1D & y){
	y.WriteIt(fout);
	return fout;
}



