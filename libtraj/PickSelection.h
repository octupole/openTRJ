/*
 * PickSelection.h
 *
 *  Created on: Jul 15, 2015
 *      Author: marchi
 */

#ifndef SRC_PICKSELECTION_H_
#define SRC_PICKSELECTION_H_
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <set>

#include <myEnums.hpp>
#include "Topol.h"

using std::map;
using std::set;
using std::vector;
using std::string;
using std::cout;
using std::cin;
using std::endl;
typedef map<string,vector<vector<int> > > mappa;

using namespace Topol_NS;
using namespace Enums;
struct ListTypeRes{
	vector<string> operator()(Topol_NS::Topol & MyTop,string myType){
		vector<string> & AllNames=MyTop.getStringTypeNames();
		vector<string> & AllRes=MyTop.getrestype_s();
		vector<string> NewRes;
		for(size_t o=0;o<AllRes.size();o++){
			if(AllRes[o].find(myType) != string::npos)
				NewRes.push_back(AllNames[o]);
		}
		return NewRes;
	}
};


struct PickSelection{
	vector<string> SelRes;
	PickSelection(vector<string> & ss): SelRes{ss}{}

	template <Enums::myAtoms enu>
	vector<string> Select(mappa & Def, Topol & MyTop){
		return __doSelection(Def,MyTop);
	}
private:
	vector<string> __doSelection(mappa & Def, Topol & MyTop){
		vector<string> CSelect;
		if(SelRes.empty()){
			cout << "\n System residues are: ";
			for(mappa::iterator it=Def.begin();it != Def.end();it++){
				cout << it->first << " " ;
			}
			cout << endl;
			cout << "\nSelection of the system residues not given. \nPick "
					"the residues now  and type ctrl-D to finish: ";
		} else {

			vector<string>::iterator it=SelRes.begin();
			vector<string> newsel;


			for(;it!=SelRes.end();it++){
				vector<string> NewRes=ListTypeRes()(MyTop,*it);
				if(NewRes.size() != 0)
					copy(NewRes.begin(),NewRes.end(),back_inserter(newsel));
				else
					newsel.push_back(*it);
			}
			SelRes=newsel;
			for(unsigned int n=0;n< SelRes.size(); n++)
				if(!Def.count(SelRes[n])) {
					if(Def.count(SelRes[n]+"_P") && Def.count(SelRes[n]+"_O")){
						auto it=SelRes.begin()+n;
						string tmp=SelRes[n]+"_O";
						SelRes[n]+="_P";
						SelRes.insert(it,tmp);
					} else
						CSelect.push_back(SelRes[n]);
				}

			if(!CSelect.empty()){
				cout << "\n System residues are: ";
				for(mappa::iterator it=Def.begin();it != Def.end();it++){
					cout << it->first << " " ;
				}
				cout << endl;
				cout << "\n Residue(s) " ;
				for(unsigned m=0; m<CSelect.size();m++)
					cout << CSelect[m] << "  ";
				cout <<" are not on the list.\n\nPick "
						"the residues now  and type ctrl-D to finish:\n";
			} else {
				return SelRes;
			}
		}
		if(SelRes.empty() || CSelect.size()){
			SelRes.clear();
			bool ok=true;
			do {
				ok=true;
				string dummy;
				vector<string> newsel;
				while(cin>>dummy) {
					vector<string> NewRes=ListTypeRes()(MyTop,dummy);
					if(NewRes.size() != 0)
						copy(NewRes.begin(),NewRes.end(),back_inserter(newsel));
					else
						newsel.push_back(dummy);
				};
				SelRes=newsel;
				vector<string> unFoundRes;
				for(unsigned int n=0;n < SelRes.size();n++){
					if(!MyTop.CheckResidue(SelRes[n])) {
						ok=false;
						unFoundRes.push_back(SelRes[n]);
					}
				}
				if(ok) {
					cout << "\nSelected Residues are:  ";
					for(unsigned int n=0;n < SelRes.size();n++)
						cout << SelRes[n] << " ";
					cout <<endl;
				} else{
					cout << "Residues \"";
					for(auto v: unFoundRes)
						cout << v << " ";
					cout << "\" not on the list. Abort.\n";
					exit(1);
				}
			} while(!ok);
			return SelRes;
		}
	return SelRes;

	}

};


#endif /* SRC_PICKSELECTION_H_ */
