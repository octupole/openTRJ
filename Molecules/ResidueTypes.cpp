/*
 * ResidueTypes.cpp
 *
 *  Created on: May 21, 2012
 *      Author: marchi
 */

#include "ResidueTypes.h"

namespace Topol_NS {
string ResidueTypes::_types[]={"Prot","Water","Ions","Surf","Dna", "Cofact","OrgSol","DetPol", "DetOil","Other"};
#include <Residues.txt>
vector<string> ResidueTypes::_notfound=vector<string>();
	ResidueTypes::~ResidueTypes() {
		// TODO Auto-generated destructor stub
	}
int ResidueTypes::find(string & y ){
	string x=y+" ";
	try{
		for(int n=0;n<TOTTYPES;n++){
			size_t found=_Residue[n].find(x.c_str());
			while(found != string::npos){
				if(found == 0) return n;
				else if(_Residue[n][found-1] == ' ') return n;
				found =_Residue[n].find(x.c_str(),found+1) ;
			}

		}
		if(y.size()>3)
			for(int n=0;n<TOTTYPES;n++){
				string tmp=y.substr(0,y.size()-1)+" ";
				if(_Residue[n].find(tmp.c_str()) != string::npos) return n;
			}
		if(std::find(_notfound.begin(),_notfound.end(),y) == _notfound.end()) {
			_notfound.push_back(y);
			throw string("\nWarning: Can\'t find pdb residue " + y + " in the internal list. Its value will be set to OTHER");
		} else return TOTTYPES-1;

	}
	catch(const string & s){
		cout << s << endl;
		return TOTTYPES-1;
	}
	return -1;
}

};


