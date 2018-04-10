/*
 * TopolPDB.cpp
 *
 *  Created on: Aug 24, 2013
 *      Author: marchi
 */

#include "TopolPDB.h"

namespace Topol_NS {

auto countRes=[](size_t & n, string & oldRes, string currRes) {
	if(oldRes == " ") {
		oldRes=currRes;
	}
	if(oldRes != currRes){
		n++;
		oldRes=currRes;
	}
};
auto myName=[](string s, string res){
	const string knownRes{"SOD NA CLA K MG CA ZN"};
	string atom;

	if(knownRes.find(s) != string::npos){
		if(s == "SOD")
			atom="Na";
		else if(s == "CLA")
			atom="Cl";
		else {
			std::transform(s.begin()+1,s.end(),s.begin()+1,tolower);
			atom=s;
		}
	} else{
		if(s.find_first_not_of("1234") == 1){
			if(s.compare(1,1,"H") == 0){
				atom="H";
			}
		}else if(s == "FE" || s == "Fe")
			atom="Fe";
		else
			atom=s.substr(0,1);
	}
	return atom;
};
void TopolPDB::operator()(const vector<string> & data0 ){
	if(data0.empty()) return;
	map<string,vector<string> > mapps;
	vector<string> data=vector<string>(data0.size(),string());
	int ia=0;
	for(size_t i=0;i<data0.size();i++){
		if(data0[i].substr(0,7).find("ATOM") == 0 || data0[i].substr(0,7).find("HETATM") == 0 ) {
			data[ia]=data0[i];
			ia++;
		}
	}
	data.resize(ia);
//	op_ResNo tmp;
//
//	std::transform(data.begin(),data.end(),data.begin(),tmp);

	size_t nresMax{0};


	string oldRes{" "};
	vector<size_t> vt;
	for(auto op: data){
		string Res=op.substr(22,5);
		countRes(nresMax,oldRes,Res);
		vt.push_back(nresMax);
	}
	vector<vector<string> > vstring=vector<vector<string> > (nresMax+1);

	for(size_t i=0;i<data.size();i++){
		size_t nres=vt[i];
		vstring[nres].push_back(data[i]);
	}
	size_t ntot=0;
	for(size_t n=0; n< vstring.size(); n++){
		ntot+=vstring[n].size();
	}

	for(size_t n=0; n< vstring.size(); n++){
		vector<PDBdata> tmp=SplitResidue(vstring[n]);
		PDB.insert(PDB.end(),tmp.begin(),tmp.end());
	}
}
void TopolPDB::sDetPolsegment(string s0,string s){
	if(s0.find("None") != string::npos) return;
	DetPolSegment={vector<string>(2)};
	DetPolSegment[0]=s0;
	vector<string> vs;
	vector<string> vs0=split(s);

	for(auto it: vs0){
		cleanString(it);
		DetPolSegment[1]+=it+" ";
	}
	string subP=s0+"_P";
	string subO=s0+"_O";
	ResidueTypes::addDetPolar(subP);
	ResidueTypes::addDetOil(subO);
}
vector<PDBdata> TopolPDB::SplitResidue(vector<string> & v0){
	string sub2=v0[0].substr(17,4); // Residue name
	sub2.erase(remove_if(sub2.begin(),sub2.end(),::isspace),sub2.end());
	vector<int> atmindx=vector<int>(v0.size());
	vector<int> resindx=vector<int>(v0.size());
	vector<string> resname=vector<string>(v0.size());

	vector<PDBdata> y=vector<PDBdata>(v0.size());
	vector<PDBdata> yout=vector<PDBdata>(v0.size());
	for(size_t o=0;o < v0.size(); o++){
		string sub1=v0[o].substr(11,5); // Atom name
		sub1.erase(remove_if(sub1.begin(),sub1.end(),::isspace),sub1.end());
		y[o].atn=sub1;
		std::stringstream ss(v0[o].substr(31,24)); // coordinates
		ss>> y[o].x;
		if(isdigit(sub1[0])) y[o].type=sub1.substr(1,4)+sub1.substr(0,1);
		else y[o].type=sub1.substr(0,5);
		// And now let us do the subunit stuff
		std::stringstream(v0[o].substr(61,6))>>y[o].domd;
		y[o].domn=v0[o].substr(21,1);
	}
	if(!DetPolSegment.empty()){
		if(sub2.find(DetPolSegment[0]) != string::npos) {
			string subP=sub2+"_P";
			string subO=sub2+"_O";
			for(size_t o=0;o < y.size(); o++){
				resname[o]=subP;
				atmindx[o]=o;
				resindx[o]=offsetRes;
				if(DetPolSegment[1].find(y[o].atn) == string::npos) {
					resname[o]=subO;
					resindx[o]=offsetRes+1;
				}
			}
			offsetRes+=2;
		} else{
			for(size_t o=0; o < v0.size();o++){
				atmindx[o]=o;
				resindx[o]=offsetRes;
				resname[o]=sub2;
			}
			offsetRes++;
		}
	} else{
		for(size_t o=0; o < v0.size();o++){
			atmindx[o]=o;
			resindx[o]=offsetRes;
			resname[o]=sub2;
		}
		offsetRes++;

	}
	try{
		for(size_t o=0;o < y.size(); o++){
			yout[o].at=atmindx[o]+offset;
			yout[o].res=resindx[o];
			yout[o].resn=resname[o];
			yout[o].atn=y[o].atn;
			yout[o].type=y[o].type;
			yout[o].x=y[o].x;
			yout[o].domd=y[o].domd;
			yout[o].domn=y[o].domn;
			yout[o].first=v0[o].substr(0,61);
			yout[o].last=v0[o].substr(66,14);
			string atm=myName(yout[o].atn,yout[o].resn);
			if(Masses.count(atm) != 1) throw string("Unknown atom found");
			yout[o].mass=Masses[atm];
		}
	}catch(const string & s){
		cout << s << endl;
		exit(1);
	}
	offset+=v0.size();

	return yout;
}

TopolPDB::~TopolPDB() {
	// TODO Auto-generated destructor stub
}

} /* namespace Topol_NS */
