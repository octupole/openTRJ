/*
 * ResidueTypes.h
 *
 *  Created on: May 21, 2012
 *      Author: marchi
 */

#ifndef RESIDUETYPES_H_
#define RESIDUETYPES_H_
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>

using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace Topol_NS {
const int TOTTYPES=10;
enum atomSpecies{Prot, Water, Ions, Surf, Dna, Cofact, OrgSol, DetPol, DetOil,Other};
struct MyResidue{
	int n;
	string name;
};
class ResidueTypes {
	static string _Residue[TOTTYPES];
	static string _types[TOTTYPES];
	static vector<string> _notfound;
	vector<string> Residue;
	vector<string> __ResTypes;
	int WaterId;
public:
	ResidueTypes(): WaterId(0){};
	static string getType(int n){return _types[n];};
	static int getTotTypes(){return TOTTYPES;}
	static void addDetPolar(string s){_Residue[DetPol]+=s+" ";}
	static void addDetOil(string s){_Residue[DetOil]+=s+" ";}
	static string & getReslist(int i){return _Residue[i];}
	static int find(string & );

	size_t Sfind(string & y){
		vector<string>::iterator p=std::find(Residue.begin(),Residue.end(),y);
		size_t n=0;
		int tm=find(y);
		if( p == Residue.end()) {
			Residue.push_back(y);
			n=Residue.size()-1;
			__ResTypes.push_back(_types[tm]);
		} else n=std::distance(Residue.begin(),p);
		if(tm == Water) WaterId=n;
		return n;
	}
	int gWaterId(){return WaterId;}
	bool testIons(string & y ){
		if(find(y) == Ions) return true;
		return false;
	}
	string strFind(const int n){

		try {
			if(size_t(n) < Residue.size()) return Residue[n];
			else throw string(" Stop here, residue number larger than expected ");
		}
		catch(const string & s){
			cout << s << endl;
			exit(-1);
		}
		return NULL;
	}
	string strFindResType(const int n){
		try {
			if(size_t(n) < __ResTypes.size()) return __ResTypes[n];
			else throw string(" Stop here, residue number larger than expected ");
		}
		catch(const string & s){
			cout << s << endl;
			exit(-1);
		}
		return NULL;
	}
	int strFindResTypeNum(const int n){
		try {
			if(size_t(n) < __ResTypes.size()) {
				for(int m=0;m < TOTTYPES; m++){
					if(_types[m].find(__ResTypes[n]) != string::npos) return m;
				}
				throw string(" Stop here, inconsistency check true ");
			}
			else throw string(" Stop here, residue number larger than expected ");
		}
		catch(const string & s){
			cout << s << endl;
			exit(-1);
		}
		return 0;
	}
	static string types(int n){return _types[n];};
	static int Size(){return TOTTYPES;};
	static int getTypes(const string & s){
		try{
			for(int n=0;n<TOTTYPES;n++){
				if(_types[n].find(s) != string::npos) return n;
			}
			throw string("Residue type not found");
		} catch(const string & b){
			cout << b << endl;
			exit(1);
		};
		return -1;
	}
	static int getTypes(const char * s){
		try{
			for(int n=0;n<TOTTYPES;n++){
				if(_types[n].find(s) != string::npos) return n;
			}
			throw "Residue type not found";
		} catch(const string & b){
			cout << b << endl;
			exit(1);
		};
		return -1;
	}
	virtual ~ResidueTypes();
	};
};
#endif /* RESIDUETYPES_H_ */
