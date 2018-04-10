/*
 * Percolation.h
 *
 *  Created on: Sep 12, 2013
 *      Author: marchi
 */

#ifndef PERCOLATION_H_
#define PERCOLATION_H_
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include "MyUtilClass.h"
#include "LCells.h"
#include <functional>

using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::string;
using std::ostream;
using std::ofstream;
using std::function;
using namespace DVECT;
using namespace MATRIX;

typedef vector<vector<int> > listcon;

//const string exclusion="C17 C18 C20 C8 C9 C11";
const string exclusion="C15 C16 C17 C18 C19 C20 C6 C7 C8 C9 C10 C11";
//const string exclusion="S OS2 OS3 OS4";
struct Comp{
	vector<string> Res;
	vector<int> No;
	friend std::ostream & operator<<(std::ostream & fout , vector<Comp> & y){

		for(size_t o=0;o<y.size();o++){
			int ntot=0;
			for(size_t p=0;p<y[o].Res.size();p++){
				ntot+=y[o].No[p];
			}
			fout << "  Cluster = "<< std::fixed << std::setw(3) << o << " Size " << ntot << " : ";
			for(size_t p=0;p<y[o].Res.size();p++){
				fout << y[o].Res[p] << "[" << y[o].No[p] << "] ";
			}
			fout << endl;
		}
		return fout;
	}
};
const double CUTOFF=1.3;
template <typename T>
class Percolation {
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
protected:
	size_t nr{0};
	vector<vector<int> > Atoms;
	vector<int> bAtoms;
	vector<int> bIndx;
	vector<string> Name;
	vector<vector<string> >pAtn;
	vector<double> pRd;
	vector<vector<int> > Contacts;
	vector<vector<int> > Clusters;

	static vector<string> Select;
	static double PercoCutoff;
	vector<string> pResn;
	vector<Comp> ClustComp;
	double Rcut{CUTOFF};
	void rCluster(size_t);
	int count{0};
	function<double(int,int)> myPercoCutoff;
	virtual void __Writeit(ostream &);
	Matrix CO,OC;

public:
	static void setPercoCutoff(double cut){PercoCutoff=cut;}
	Percolation(){};
	Percolation(const listcon & , const vector<string> &, const vector<double> &);
	Percolation(const vector<vector<int> > & , const vector<double> &, const vector<string> &);
	Percolation(const vector<vector<int> > & , const vector<double> &, const  vector<string> &, const vector<string> &);
	void setRcut(double y){Rcut=y;};
	void Accumulate(){count++;}
	void doContacts(vector<DDvect<T>> &);
	void doContacts(vector<DDvect<T>> &, MMatrix<T> &, MMatrix<T> &);
	int gCluster();
	const listcon & getCluster() const {return Clusters;}
	const vector<int> & getFirstCluster() const {return Clusters[0];}
	listcon & getAtoms(){return Atoms;}
	listcon & getContacts(){return Contacts;}
	vector<Comp> & getClustComp(){return ClustComp;}
	vector<string> & getpResn(){return pResn;}
	virtual ~Percolation();
	friend ostream & operator<<(ostream & fout, Percolation & perco){
		perco.__Writeit(fout);
		return fout;
	}

};

#endif /* PERCOLATION_H_ */
