/*
 * TopolPDB.h
 *
 *  Created on: Aug 24, 2013
 *      Author: marchi
 */

#ifndef TOPOLPDB_H_
#define TOPOLPDB_H_

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cctype>
#include "MyUtilClass.h"
#include "Split.h"
#include "ResidueTypes.h"

namespace Topol_NS {
using namespace DVECT;

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::map;
using std::stringstream;
using std::copy;
using Dvect=DDvect<double>;

struct PDBdata{
	int at, res;
	string atn,resn,domn,type;
	Dvect x;
	float domd;
	double mass;
	string first, last;
};

typedef map<string,map<string,vector<string> > > Mymapmap;
typedef map<string,vector<string > > Mymap;
class TopolPDB {
	struct op_ResNo{
		int n;
		int oldno;
		op_ResNo(): n(0), oldno(-1){};
		string operator()(const string & y){
			string out=y;
			int currno;
			std::stringstream(out.substr(22,4))>> currno; // Residue number
			if(oldno != currno) n++;
			oldno=currno;
			std::stringstream ss;ss<<std::fixed << std::setw(5) << n;
			out.replace(22,4,ss.str());
			return out;
		}
		int getSize() const {return oldno;};
	};
	map<string,double> Masses{
			{"H",     1.00784}
			,{"He",     4.00260}
			,{"Li",     6.93800}
			,{"Be",     9.01218}
			,{"B",    10.80600}
			,{"C",    12.00960}
			,{"N",    14.00643}
			,{"O",    15.99903}
			,{"F",    18.99840}
			,{"Ne",    20.17970}
			,{"Na",    22.98977}
			,{"Mg",    24.30400}
			,{"Al",    26.98154}
			,{"Si",    28.08400}
			,{"P",    30.97376}
			,{"S",    32.05900}
			,{"Cl",    35.44600}
			,{"Ar",    39.94800}
			,{"K",    39.09830}
			,{"Ca",    40.07800}
			,{"Sc",    44.95591}
			,{"Ti",    47.86700}
			,{"V",    50.94150}
			,{"Cr",    51.99610}
			,{"Mn",    54.93804}
			,{"Fe",    55.84500}
			,{"Co",    58.93319}
			,{"Ni",    58.69340}
			,{"Cu",    63.54600}
			,{"Zn",    65.38000}
			,{"Ga",    69.72300}
			,{"Ge",    72.63000}
			,{"As",    74.92159}
			,{"Se",    78.97100}
			,{"Br",    79.90100}
			,{"Kr",    83.79800}
			,{"Rb",    85.46780}
			,{"Sr",    87.62000}
			,{"Y",    88.90584}
			,{"Zr",    91.22400}
			,{"Nb",    92.90637}
			,{"Mo",    95.95000}
			,{"Ru",   101.07000}
			,{"Rh",   102.90550}
			,{"Pd",   106.42000}
			,{"Ag",   107.86820}
			,{"Cd",   112.41400}
			,{"In",   114.81800}
			,{"Sn",   118.71000}
			,{"Sb",   121.76000}
			,{"Te",   127.60000}
			,{"I",   126.90447}
			,{"Xe",   131.29300}
			,{"Cs",   132.90545}
			,{"Ba",   137.32700}
			,{"La",   138.90547}
			,{"Ce",   140.11600}
			,{"Pr",   140.90766}
			,{"Nd",   144.24200}
			,{"Sm",   150.36000}
			,{"Eu",   151.96400}
			,{"Gd",   157.25000}
			,{"Tb",   158.92535}
			,{"Dy",   162.50000}
			,{"Ho",   164.93033}
			,{"Er",   167.25900}
			,{"Tm",   168.93422}
			,{"Yb",   173.04500}
			,{"Lu",   174.96680}
			,{"Hf",   178.49000}
			,{"Ta",   180.94788}
			,{"W",   183.84000}
			,{"Re",   186.20700}
			,{"Os",   190.23000}
			,{"Ir",   192.21700}
			,{"Pt",   195.08400}
			,{"Au",   196.96657}
			,{"Hg",   200.59200}
			,{"Tl",   204.38200}
			,{"Pb",   207.20000}
			,{"Bi",   208.98040}
	};

	vector<PDBdata> PDB;
	vector<PDBdata> SplitResidue(vector<string> &);
	int offset;
	int offsetRes;
	vector<string> DetPolSegment;
public:
	TopolPDB(): offset(0), offsetRes(0){};
	void operator()(const vector<string> &);
	void sDetPolsegment(string,string);
	PDBdata & operator[](size_t n){return PDB[n];}
	size_t Size() const {return PDB.size();}
	virtual ~TopolPDB();
};

} /* namespace Topol_NS */
#endif /* TOPOLPDB_H_ */
