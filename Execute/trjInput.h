/*
 * trjInput.h
 *
 *  Created on: May 22, 2012
 *      Author: marchi
 */

#ifndef TRJINPUT_H_
#define TRJINPUT_H_
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "../Execute/ClearUsage.h"
using std::map;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::set;

namespace trj {
enum myinput {o,dcd,xtc,b,e,skip,pdb,pVol,noHyd,Hyd,nodel,del,rd,nord,fab,nofab,test,help};
class trjInput {
protected:
	map<string,vector<string> > inmap;
	vector<string> unknownC;
	map<int,string> Usage;
public:
	trjInput(int nv,char ** v,ClearUsage &);
	virtual ~trjInput();
	vector<string> & operator[](const char * y){return inmap[y];}
	vector<string> & bTest(){return unknownC;};
	vector<string> getUsage();
};

} /* namespace trj */
#endif /* TRJINPUT_H_ */
