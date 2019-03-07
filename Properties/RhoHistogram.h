/*
 * RhoHistogram.h
 *
 *  Created on: Mar 7, 2019
 *      Author: marchi
 */

#ifndef PROPERTIES_RHOHISTOGRAM_H_
#define PROPERTIES_RHOHISTOGRAM_H_
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include "Units.hh"

using namespace Units;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::ostream;

namespace Properties {

class RhoHistogram {
	static double rcut;
	static double dx;
	vector<double> Histo;
	vector<double> Avg;
	vector<size_t> nHisto;
	double Rc{-1},mass{1.0};
	string header;
	static string fileout;
	void cPrint(ostream &);
	RhoHistogram(){};
	size_t nstep{0};
	double conv_Angstroms{10.0};
public:
	RhoHistogram(double);
	static void SetRcut(double x){rcut=x;}
	static void SetDx(double x){dx=x;}
	static void setFilename(string s){fileout=s;}

	static double getRcut(){return rcut;}
	static double getDx(){return dx;}
	static string getFilename(){return fileout;}
	void setMass(double m){mass=m;}
	void operator()(double radius, double value);
	vector<double> operator[](size_t);
	void operator()(double radius);
	RhoHistogram & operator++();
	RhoHistogram operator++(int);
	void Averages();
	void Averages(double);
	friend ostream & operator << (ostream & fout , RhoHistogram  & y ){
		y.cPrint(fout);
		return fout;
	};
	virtual ~RhoHistogram();
};

} /* namespace Properties */

#endif /* PROPERTIES_RHOHISTOGRAM_H_ */
