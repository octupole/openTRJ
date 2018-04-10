/*
 * histograms.hpp
 *
 *  Created on: Aug 16, 2013
 *      Author: marchi
 */

#ifndef HISTOGRAMS_HPP_
#define HISTOGRAMS_HPP_
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

#include <Ftypedefs.h>



using namespace std;

class hist1D{
protected:
public:
	double hisX{0};
	double nhis{0.0};
	hist1D():hisX(0.0),nhis(0){};
	hist1D(double his0,double nhis0): hisX(his0),nhis(nhis0){};
	hist1D(double his0): hisX{his0}{};
	hist1D & operator+=(hist1D a){
		hisX+=a.hisX;
		return *this;
	}
	hist1D & operator/=(double a){
		hisX/=a;
		return *this;
	}
	hist1D & operator*=(double a){
		hisX*=a;
		return *this;
	}
	double gHisx(){return hisX;}
	bool isZero(){
		return nhis == 0;
	}
	double Ratio(){
		try{
			if(isZero()) throw string("Non accumulated histogram. Something went wrong.");
			return hisX/nhis;
		}catch(const string & s){
			cout << s <<endl;
			exit(1);
		}
	}
	double gHis(){return hisX;}
	double gNhis(){return nhis;}
	hist1D & operator+=(double n){
		nhis+=n;
		return *this;
	}
	hist1D & operator++(){
		nhis++;
		return *this;
	}
	hist1D operator++(int){
		hist1D temp{*this};
		++*this;
		return temp;
	}
	hist1D & operator=(const hist1D & y){
		hisX=y.hisX;
		nhis=y.nhis;
		return *this;
	}
};
class Histogram1D{
protected:
	virtual void WriteIt(ostream & );
	vector<hist1D> hist;
public:
	double dx{0.0},cut{0.0};
	size_t HisX{0};
	string Label{" "};
	Histogram1D(){};
	Histogram1D(double,double);
	size_t Size(){return HisX;}
	virtual void operator()(double,double);

	virtual void setLabel(string s){Label=s;}

	virtual hist1D & operator[](size_t);

	Histogram1D & operator*=(double );
	Histogram1D & operator/=(double );
	Histogram1D & operator*=(size_t );
	Histogram1D & operator=(const Histogram1D &);
	bool operator==(const Histogram1D &);

	virtual Histogram1D & operator+=(const Histogram1D &);
	Histogram1D & operator+=(double n){for(size_t o=0;o<hist.size();o++){hist[o]+=n;};return *this;};
	Histogram1D & operator++(){for(size_t o=0;o<hist.size();o++){hist[o]++;};return *this;};
	Histogram1D operator++(int){Histogram1D temp{*this};++*this;return temp;};
	virtual void clear();
	friend ostream & operator<<(ostream &, Histogram1D & );
	virtual ~Histogram1D(){};
};
class PairCorr1D: public Histogram1D{
	size_t nMol{0};
	virtual void WriteIt(ostream &);
public:
	PairCorr1D(){}
	PairCorr1D(double x,double y,size_t n): Histogram1D(x,y), nMol{n}{};
	void operator()(double x,double y,size_t n){
		Histogram1D::operator ()(x,y);
		nMol=n;
	}

	friend ostream & operator<<(ostream &, PairCorr1D & );
};
#endif /* HISTOGRAMS_HPP_ */
