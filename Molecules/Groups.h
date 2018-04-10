/*
 * Groups.h
 *
 *  Created on: Jun 5, 2013
 *      Author: marchi
 */

#ifndef GROUPS_H_
#define GROUPS_H_
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include "MyUtilClass.h"
#include "ResidueTypes.h"
using namespace DVECT;
using Dvect=DDvect<double>;
using std::vector;
namespace Topol_NS {
struct MyData {
	Dvect v;
	string name;
	MyData(Dvect & v0, string & name0):v(v0), name(name0) {};
};
class Groups {
	vector<bool> MyBeta;

public:
	Groups();
	Groups(vector<Dvect> &,vector<MyResidue> &, vector<string> &);
	void Init(vector<Dvect> &,vector<MyResidue> &, vector<string> &);
	void operator()(vector<Dvect> &,vector<MyResidue> &, vector<string> &);
	bool getBeta(size_t i){
		try{
			if(i > MyBeta.size()) throw string(" In groups MyBeta smaller than expected ");
		} catch(const string & s) {
			cout << s << endl;
		}
		return MyBeta[i];
	};
	virtual ~Groups();
};
} /* namespace Topol_NS */
#endif /* GROUPS_H_ */
