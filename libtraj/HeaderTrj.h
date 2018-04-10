/*
 * HeaderTrj.h
 *
 *  Created on: May 13, 2012
 *      Author: marchi
 */

#ifndef HEADERTRJ_H_
#define HEADERTRJ_H_
#include <vector>
#include <string>

#include "../libtraj/FstreamC.h"
#include "../libtraj/FstreamF.h"

using std::vector;
using std::string;
using std::ifstream;
using std::ios;

const int HDIM=4;
const int TDIM=80;
const int FORTRANBYTES=4;
class HeaderTrj {
	char hdr[HDIM];
	int nfr{0},istart{0},natoms{0};
	vector<string> title;
	void ReadHeader(FstreamC *);
	void ReadHeader(FstreamF *);
	void ReadHeader(std::ifstream &);
public:
	HeaderTrj();
	virtual ~HeaderTrj();
	int getNFR(){return nfr;}
	int getNatoms(){return natoms;}
	bool check(int n){return n==natoms;}
	friend Fstream & operator>>(Fstream & fin, HeaderTrj & y){
		if(FstreamC * finC=dynamic_cast<FstreamC *> (&fin))
			y.ReadHeader(finC);
		else if(FstreamF * finF=dynamic_cast<FstreamF *> (&fin))
			y.ReadHeader(finF);
		return fin;
	}

	friend std::ifstream & operator>>(std::ifstream & fin, HeaderTrj & y){
		y.ReadHeader(fin);
		return fin;
	}
};

#endif /* HEADERTRJ_H_ */
