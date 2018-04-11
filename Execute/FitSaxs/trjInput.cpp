/*
 * trjInput.cpp
 *
 *  Created on: May 22, 2012
 *      Author: marchi
 */

#include "trjInput.h"

namespace trj {
trjInput::trjInput(int ntot,char ** v) {
	auto ClrU=[this](std::initializer_list<int> list){for(auto op: list) Usage[op].clear();};
	vector<string> in;

	inmap["-o"]=in;
	inmap["-in"]=in;

	map<string,vector<string> >::iterator it=inmap.begin();
	for(int n=0;it!=inmap.end();++it,n++){
		Usage[n]=" ";
	}
	Usage[0]="\t -o fileout: The new computed function fitted to experiment\n";
	Usage[1]="\t -in Read the computed and experimental function \n";

	Usage[2]="\t -help // write some on line help \n";

	vector<string> vv0;
	string vv;
	for(auto o=1;o<ntot;o++) vv0.push_back(v[o]);

	string key;
	for(auto tmp0: vv0){
		if(tmp0[0] =='-'){
			key.assign(tmp0);
			if(inmap.find(key) != inmap.end()){
				if(inmap[key].empty()) inmap[key].push_back(tmp0);
			} else{
				unknownC.push_back(key);
				inmap[key].push_back(tmp0);
			}
		}
		else{
			inmap[key].push_back(tmp0);
		}
	}
}

vector<string> trjInput::getUsage(){
	vector<string> use;
	auto it=Usage.begin();
	for(;it!=Usage.end();++it){
		if(!Usage.empty())
			use.push_back(it->second);
	}
	return use;
}
trjInput::~trjInput() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */
