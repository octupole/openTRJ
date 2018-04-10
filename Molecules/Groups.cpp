/*
 * Groups.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: marchi
 */

#include "Groups.h"
namespace Topol_NS {
bool isPolar(string & a){if(a.substr(0,1) == "N" || a.substr(0,1) == "O") return true; return false;};
bool isH(string & a){if(a.substr(0,1) == "H") return true; return false;};
bool isC(string & a){if(a.substr(0,1) == "C") return true; return false;};

Groups::Groups() {
	// TODO Auto-generated constructor stub

}
Groups::Groups(vector<Dvect> & xc,vector<MyResidue> & res, vector<string> & name){
	Groups::Init(xc,res,name);
}
void Groups::operator()(vector<Dvect> & xc,vector<MyResidue> & res, vector<string> & name){
	Groups::Init(xc,res,name);
}
void Groups::Init(vector<Dvect> & xc,vector<MyResidue> & res, vector<string> & name){
	int natoms=res.size();
	vector<int> Tag(natoms);
	for(int i=0;i<natoms;i++){
		if(i==0) Tag[i]=0;
		else {
			Tag[i]=Tag[i-1];
			if(res[i-1].n != res[i].n) Tag[i]=Tag[i-1]+1;
		}
	}
	int resDim=*max_element(Tag.begin(),Tag.end())+1;
	vector<vector<MyData> > MyD(resDim);
	for(int i=0;i<natoms;i++){
		MyD[Tag[i]].push_back(MyData(xc[i],name[i]));

	}

	for(int n=0;n<resDim;n++){
		vector<MyData> atm=MyD[n];
		int natom=atm.size();
		vector<vector<int> > nei(natom);
		for(int i=0;i<natom;i++){
			Dvect xi=atm[i].v;
			for(int j=0;j<natom;j++){
				if(i!=j){
					double dist=xi.Dist(atm[j].v);
					if(dist < 1.6) {
						nei[i].push_back(j);
					}
				}
			}
		}
		vector<bool> isPhob(natom);
		for(int i=0;i<natom;i++)
			isPhob[i]=false;
		for(int i=0;i<natom;i++){
			string element_i=atm[i].name.substr(0,1);
			string name_i=atm[i].name;
			if(element_i == "C" && name_i != "CA"){
				vector<string> name_int;
				for(size_t j=0;j<nei[i].size();j++){
					string name_j=atm[nei[i][j]].name;
					name_int.push_back(name_j);
				}
				int NoH=count_if(name_int.begin(),name_int.end(),isH);
				int NoP=count_if(name_int.begin(),name_int.end(),isPolar);
				if(! NoP){
					isPhob[i]=true;
					if (NoH){
						for(size_t j=0;j<nei[i].size();j++){
							if(atm[nei[i][j]].name.substr(0,1) == "H") isPhob[nei[i][j]]=true;
						}
					}
				}
			}
		}
		for(int i=0;i<natom;i++){
			MyBeta.push_back(isPhob[i]);
		}
	}
}

Groups::~Groups() {
	// TODO Auto-generated destructor stub
}
} /* namespace Topol_NS */
