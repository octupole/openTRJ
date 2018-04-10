/*
 * SaxsPadding.h
 *
 *  Created on: Feb 26, 2016
 *      Author: marchi
 */

#ifndef SRC_SAXSPADDING_H_
#define SRC_SAXSPADDING_H_
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iomanip>
#include "Finalize.h"



using namespace std;
/** \brief The class selects what atom to use in the padding of the supercell
 *
 */
class SaxsPadding {
	vector<string> List;
	vector<double> iList;
	map<string,double> MapResidue;
	bool Have_CheckedIt{false};
public:
	SaxsPadding(){};
	SaxsPadding(vector<string> & vstr,vector<double> vint){
		this->operator()(vstr,vint);
	};
	bool Have_MyPadding(){if(List.size()) {return true;}return false;}
	vector<string> & gList(){return List;}
	vector<double> & giList(){return iList;}
	void operator()(vector<string> & vstr,vector<double> vint);
	void setMapResidue(const map<string,vector<string> > &);
	const map<string,double> & getMapResidue() const {return MapResidue;}
	void CheckIt(vector<string> & Ref, vector<string> & Solu);
	virtual ~SaxsPadding();
};

#endif /* SRC_SAXSPADDING_H_ */
