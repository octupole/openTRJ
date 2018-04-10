/*
 * Topol.h
 *
 *  Created on: May 12, 2012
 *      Author: marchi
 */

#ifndef TOPOL_H_
#define TOPOL_H_
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sstream>
#include <map>
#include <numeric>
#include "ResidueTypes.h"
#include "MyUtilClass.h"
#include "Groups.h"
#include "TopolPDB.h"
#include "myEnums.hpp"
//#include "Split.h"
#include "Finalize.h"

using std::accumulate;

using namespace DVECT;

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::map;
using std::pair;
using Dvect=DVECT::DDvect<double>;
using namespace Enums;
namespace Topol_NS {
enum {Slt, Sol};

class MySubunit{
	bool b;
	string c;
	int n;
public:
	MySubunit(){b=false;c=' ';n=-1;};
	MySubunit(bool b1, string c1): n(-1) {
		b=b1;c=c1;
	}
	bool isHydrophobic(MySubunit a){
		if(c != a.c &&(b==true && a.b == true)) return true;
		return false;
	}
	bool isHydrophilic(MySubunit a){
		if(c != a.c &&(b==false && a.b == false)) return true;
		return false;
	}
	bool isMixed(MySubunit a){
		if(c != a.c) {
			if(b == true && a.b == false) return true;
			if(b == false && a.b == true) return true;
		}
		return false;
	}

	void putN(int m){n=m;};
	void putB(bool b0){b=b0;};
	void putC(string c0){c=c0;};
	string getC(){return c;};
	bool getB(){return b;};
	int getN(){return n;};
};

class Topol {
protected:
	size_t nr{0};
	size_t nres{0};
	size_t maxtype{0};
	int bGMax{0};

	vector<vector<int> > RefRes,SelRes;
	vector<string> atomname;
	vector<int> atomtypeNo;
	vector<string> resinfo;
	vector<int> resind;
	vector<vector<int> > cidx;
	vector<string> TypeNames;
	vector<double> rd;
	vector<double> Mass;
	vector<string> atSFactor;
	vector<int> bG;
	vector<int> ResType;
	vector<int> atResType;
	vector<string> restype_s;
	vector<int> restype_i;
	vector<int> ReferenceResidues;
	map<int,vector<int> > MCIndex;
	map<string,vector<vector<int> > > DefRes;
	map<string,vector<vector<int> > > DefResNH;

	vector<vector<int> > vG;
	vector<MySubunit> MySub;
	map<string,vector<int> > DefDomain;
	vector<vector<int> > CIndex;
	map<string,vector<string> > MapElements; // Map of element per residue

	void ExtractInfo(TopolPDB & ,bool);
	void ExtractInfoMic(TopolPDB & ,bool);

	void FindGroups(vector<Dvect> &,int);
	void AllignGroups(vector<Dvect> &,int &,vector<int>&,int );
	TopolPDB PDB;
public:
	static std::ofstream * fDomainPDB;
	static bool bAuto;
	static bool bDefDomain;
	Topol();
	Topol(TopolPDB &, bool);
	Topol & operator()(TopolPDB &,bool);
	virtual ~Topol();
	bool CheckResidue(const string & s);
	vector<string> & getStringTypeNames(){return TypeNames;};
	TopolPDB & gPDB(){return PDB;};
	map<string,vector<vector<int> > > & getDef(){return DefRes;}
	map<string,vector<vector<int> > > & getDefNH(){return DefResNH;}
	vector<int> & gReferenceResidues(){return ReferenceResidues;};
	void sReferenceResidues(vector<int> & v){ReferenceResidues=v;};
	void InitSelection(vector<string> &, Enums::myAtoms);

	size_t Size(){return nr;}
	size_t ResSize(){return nres;}
	int ResIndex(const int n){return resind[n];};
	vector<int> & getResind(){return resind;};
	vector<int> & gatResType(){return atResType;}
	unsigned int getCidxSize(){return cidx.size();}
	const vector<int> & getCidx(const size_t n){return cidx[n];}

	const vector<string> & getAtomName(){return atomname;}
	const vector<int> & getAtomTypeNo(){return atomtypeNo;}
	vector<string> getTypeNames(){return TypeNames;};
	void setrdtoone(){for(unsigned int n=0;n<rd.size();n++) rd[n]=1.0;}
	vector<string> & getrestype_s(){return restype_s;};
	vector<int> & getResType(){return ResType;};
	vector<int> & getrestype_i(){return restype_i;};
	vector<vector<int> > & gCIndex(){return CIndex;}
	string AtomInfo(size_t n){
		try {if(n>=nr) throw "Atom number beyond bounds ";}
		catch(const char * s){cout << s << n << " larger than " << nr << endl;exit(1);}
		return atomname[n];
	}
	vector<string> & getResinfo(){return resinfo;};
	string AtomResidue(int n){
		try {if(static_cast<size_t> (n) >=nr) throw "Atom number beyond bounds ";}
		catch(const char * s){cout << s << n << " larger than " << nr << endl;exit(1);}
		return resinfo[resind[n]];
	}
	string ResInfo(size_t n){
		try {if(n>=nres) throw "Residue number beyond bounds ";}
		catch(const char * s){cout << s << endl;exit(1);}
		return resinfo[n];
	}
	const vector<double> & getrd() const {return rd;};
	const vector<double> & getMass() const {return Mass;};
	const vector<string> & get_atSFactor() const {return atSFactor;}
	vector<MySubunit> & getMySub(){
		return MySub;
	};
	const map<string,vector<string> > & gMapElements() const {return MapElements;}
	int getmaxt(){return maxtype;}
	const vector<int> & gbG(){return bG;};
	int gbGMax(){return bGMax;}
	const vector<vector<int> > & gvG(){return vG;};
	const vector<int> & nAtomRes(size_t o){return MCIndex[o];};
	void AddDomain(const string & a, const vector<int> & b){
		DefDomain[a]=b;
	}
};

} /* namespace Topol_NS */
#endif /* TOPOL_H_ */
