/*
 * Topol.cpp
 *
 *  Created on: May 12, 2012
 *      Author: marchi
 */

#include "Topol.h"

namespace Topol_NS {
bool myunique(MyResidue & one, MyResidue & two){
	return (one.n==two.n);
}

std::ofstream * Topol::fDomainPDB=NULL;
bool Topol::bAuto=false;
bool Topol::bDefDomain=false;

struct op_sort{
	bool operator()(MyResidue x, MyResidue y){return (x.n < y.n);};
};
class my_radius{
	const int nn;
	enum {H,C,N,O,S,P,Na,Cl,Mg,Ca};
	vector<double> radius;
	vector<string> rtype;
public:
	my_radius():nn(10){
		double tmp[]={0.03,0.175,0.165,0.155,0.176,0.19,0.12,0.202,0.105,0.12};
		string tmp2[]={"H","C","N","O","S","P","NA","CL","MG","CA"};
		radius=vector<double>(tmp,tmp+nn);
		rtype=vector<string>(tmp2,tmp2+nn);
	}
	double operator()(const string & x){
		try{
			if(x.find_first_not_of("1234") == 1){
				if(x.compare(1,1,"H") == 0){
					return radius[0];
				}
			} else
				if(x.compare(0,1,"H") == 0){
					return radius[0];
				}

			for(size_t n=0;n<rtype.size();n++)
				if(x == rtype[n]) return radius[n];


			for(size_t n=0;n<rtype.size();n++)
				if(x.substr(0,1) ==rtype[n]) return radius[n];
			throw string("Stop Here: Can\'t find radius to assign to element type "+ x);
		}
		catch(const string & s){
			cout << s<< endl;
			Finale::Finalize::Final();
			exit(-1);
		}
		return -1.0;
	}
};
class op_compare{
	int seq, indx, offset;
	string sindx;
public:
	op_compare(int nseq,int nindx):seq(nseq), indx(nindx), offset(0),sindx(" ") {};
	op_compare():seq(0),indx(-1),offset(0), sindx(" "){};
	MyResidue operator()(MyResidue m){
		int n=m.n+offset;
		string ns=m.name;
		if(indx==-1){
			indx=n;
			sindx=ns;
		}
		else if(indx != n){
			seq++;
			indx=n;
			sindx=ns;
		} else if(sindx != ns){
			seq++;
			offset++;
			sindx=ns;
			indx=n+1;
		}
		m.n=seq;
		m.name=sindx;
		return m;
	}
};

class op_index{
	int seq, indx;
public:
	op_index(int nseq,int nindx):seq(nseq), indx(nindx){};
	op_index():seq(1),indx(-1){};
	int operator()(int m){
		if(!m) return 0;
		int n=m;
		if(indx==-1)indx=n;
		else if(indx != n){
			seq++;
			indx=n;
		}
		return seq;
	}
};
Topol::Topol(){}
Topol::Topol(TopolPDB & y, bool x){
	ExtractInfo(y,x);
	ExtractInfoMic(y,x);
}
Topol & Topol::operator()(TopolPDB & y,bool x){
	ExtractInfo(y,x);
	ExtractInfoMic(y,x);

	return *this;
}

Topol::~Topol() {
	// TODO Auto-generated destructor stub
}
bool Topol::CheckResidue(const string & y){
	return DefRes.count(y) > 0;
}
void Topol::ExtractInfoMic(TopolPDB & data,bool bRD){
	int nat=0;
	MyResidue mytmp;
	map<string,map<int,vector<int> > > Defs,DefsNH;

	for(size_t i=0;i<data.Size();i++){
		string sub1=data[i].atn; // Atom name
		string sub2=data[i].resn; // Residue name
		int nn=data[i].res;
		Defs[sub2][nn].push_back(nat);
		if(sub1.find_first_not_of("1234") == 1){
			if(sub1.compare(1,1,"H") != 0){
				DefsNH[sub2][nn].push_back(nat);
			}
		} else
			if(sub1.compare(0,1,"H") != 0){
				DefsNH[sub2][nn].push_back(nat);
			}

		nat++;
	}
	typedef map<string,map<int,vector<int> > > mappa;
	mappa::iterator itt=Defs.begin();
	for(;itt != Defs.end();itt++){
		map<int,vector<int> > & idx=itt->second;
		DefRes[itt->first]=vector<vector<int> >(idx.size());
		map<int,vector<int> >::iterator ipp=idx.begin();
		int ia=0;
		for(;ipp != idx.end() ;ipp++){

			DefRes[itt->first][ia]=ipp->second;
			ia++;
		}
	}
	itt=DefsNH.begin();
	for(;itt != DefsNH.end();itt++){

		map<int,vector<int> > & idx=itt->second;
		DefResNH[itt->first]=vector<vector<int> >(idx.size());
		map<int,vector<int> >::iterator ipp=idx.begin();
		int ia=0;
		for(;ipp != idx.end() ;ipp++){
			DefResNH[itt->first][ia]=ipp->second;
			ia++;
		}
	}
}

void Topol::AllignGroups(vector<Dvect> & y, int & ng,vector<int>& bG,int tname){
	const double CUT0=1.3;
	const double CUT1=1.5;
	for(unsigned int m=0; m<cidx[tname].size(); m++){
		int n=cidx[tname][m] ;
		if(!atomname[n].compare(0,1,"HXXXXX",0,1) || bG[n]) continue;
		ng++;
		for(unsigned int p=0; p<cidx[tname].size(); p++){
			int o=cidx[tname][p] ;
			if(atomname[o].compare(0,1,"HXXXXX",0,1) || n == o || bG[o]) continue;
			if(y[n].Dist(y[o]) < CUT0){
				bG[n]=ng;
				bG[o]=ng;
			}
		}
	}
	for(unsigned int m=0; m<cidx[tname].size(); m++){
		int n=cidx[tname][m] ;
		if(bG[n]) continue;
		ng++;
		for(unsigned int p=0; p<cidx[tname].size(); p++){
			int o=cidx[tname][p] ;
			if(bG[o]|| n == o) continue;
			if(y[n].Dist(y[o]) < CUT1 ){
				bG[n]=ng;
				bG[o]=ng;
			}
		}
		if(bG[n]) continue;
		bG[n]=ng;
	}
}
void Topol::FindGroups(vector<Dvect> & y, int Water){
	int ng=0;
	bG.clear();
	bG=vector<int>(y.size(),0);
	for(unsigned int n=0;n<cidx.size();n++){
		if(int(n) != Water) AllignGroups(y,ng,bG,n);
	}
 	std::transform(bG.begin(),bG.end(),bG.begin(),op_index());
 	for(unsigned int m=0;m< cidx[Water].size();m++) bG[cidx[Water][m]]=0;
 	bGMax=*max_element(bG.begin(),bG.end());
 	vG=vector<vector<int> >(bGMax+1);
	for(unsigned int n=0;n<bG.size();n++)
		vG[bG[n]].push_back(n);
}

void Topol::InitSelection(vector<string> & y, Enums::myAtoms MM)
{
	vector<vector<int> > * mySel;

	try{
		switch(MM){
		case Selection:
			mySel=&SelRes;
			break;
		case Reference:
			mySel=&RefRes;
			break;
		default:
			throw string(" ");
		}
	}catch(const string & s){cout << s <<endl;Finale::Finalize::Final();exit(1);}

	vector<int> Type=atomtypeNo;
	try{
		for(size_t o=0;o<y.size();o++){
			if(!CheckResidue(y[o])) throw string(" Selected residue " + y[o] + " does not exist in the system. Abort");
		}
	} catch(const string & s){
		cout << s << endl;
		Finale::Finalize::Final();
		exit(1);
	}
	string str=std::accumulate(y.begin(),y.end(),string(""),[](string o, string t){return o+" "+t;});
	vector<int> tmp;
	for(size_t o=0;o< nres;o++){
		if(str.find(resinfo[o]) != string::npos){
			mySel->push_back(MCIndex[o]);
			if(MM == Reference){
				tmp.push_back(o);
			}
		}
	}
	ReferenceResidues=tmp;
}

void Topol::ExtractInfo(TopolPDB & data,bool bRD){
	atomname.clear();
	resinfo.clear();
	MySub.clear();
	CIndex.clear();
	PDB=data;

	MyResidue mytmp;
	vector<MyResidue> atres;
	vector<Dvect> xc;
	vector<int> resn;
	ResidueTypes mytype;
	map<int,vector<int> > cidx0;
	map<string,int> MapRes;
	vector<string> domain;
	int nat=0;

	for(size_t i=0;i<data.Size();i++){
		string sub1=data[i].atn; // Atom name
		string sub2=data[i].resn; // Residue name
		string type=data[i].type; // Atom type

		mytmp.n=data[i].res;
		if(MapRes.find(sub2) == MapRes.end()){
			MapRes[sub2]=mytmp.n;
		}

		xc.push_back(data[i].x);
		Mass.push_back(data[i].mass);
		int ntype=mytype.Sfind(sub2);
		cidx0[ntype].push_back(nat);
		atomtypeNo.push_back(mytype.strFindResTypeNum(ntype));
		MCIndex[data[i].res].push_back(nat);
// And now let us do the subunit stuff
		float Mya=data[i].domd;
		string Myc=data[i].domn;
		MySubunit MySub0(Mya!=0.0, Myc);
		MySub.push_back(MySub0);
		domain.push_back(Myc);
// Sub unit stuff
		string sub4;
		if(::isdigit(sub1[0])) {
			sub1=sub1.substr(1,sub1.size()-1)+sub1.substr(0,1);
		}
		sub4=sub1.substr(0,2);
		mytmp.name=sub2;
		atomname.push_back(sub1);
		double myrd=0.193;
		if(bRD) rd.push_back(my_radius()(sub4));
		else rd.push_back(myrd);
		string type1{0};
		if(mytype.strFindResType(ntype) == "Prot"){
			type1=type.substr(0,1);
			if(mytype.strFind(ntype).find("HEM") != std::string::npos && sub1.find("FE") != string::npos){
				string tmp=type.substr(1,type.size()-1);
				std::transform(tmp.begin(),tmp.end(),tmp.begin(),::tolower);
				type1=type.substr(0,1)+tmp;
			}

		} else if(mytype.strFindResType(ntype) == "Ions"){
			string tmp=type.substr(1,type.size()-1);
			std::transform(tmp.begin(),tmp.end(),tmp.begin(),::tolower);
			type1=type.substr(0,1)+tmp;
		} else{
			type1=type.substr(0,1);
		}

		atSFactor.push_back(type1);
		atres.push_back(mytmp);
		nat++;
	}

	for(auto it=MapRes.begin(); it != MapRes.end();++it){
		for(size_t o=0;o<MCIndex[it->second].size();o++){
			auto o1=atSFactor[MCIndex[it->second][o]];
			MapElements[it->first].push_back(o1);
		}
	}

	auto im=MCIndex.begin();
	nat=0;
	CIndex=vector<vector<int> >(MCIndex.size());
	for(;im != MCIndex.end();im++){
		CIndex[nat++]=im->second;
	}

	if(bDefDomain) {
		map<string,vector<int> >::iterator it=DefDomain.begin();
		vector<bool> test0(MySub.size());
		std::fill(test0.begin(),test0.end(),false);
		for(;it!=DefDomain.end();++it){
			for(int n=it->second[0];n<=it->second[1];n++){
				try{
					if(!test0[n]) {
						MySub[n].putC(it->first);
						test0[n]=true;
					}
					else throw string("Domains are overlapping. Change definition.");
				}
				catch(const string & s) {
					cout << s << endl;
					Finale::Finalize::Final();
					exit(1);
				}
			}
		}
		try{if(std::count(test0.begin(),test0.end(),false)) throw string("Incomplete Domain definition! ");}
		catch(const string & s){ cout << s << endl;Finale::Finalize::Final();exit(1);}
	}
 // Interdomain stuff

 	Groups MyGroup;
 	if(bAuto) MyGroup(xc,atres,atomname);

 	vector<string>::iterator it2;
 	it2=unique(domain.begin(),domain.end());
 	domain.resize(distance(domain.begin(),it2));
 	vector<string> domain2=domain;
 	sort(domain2.begin(),domain2.end());
 	it2=unique(domain2.begin(),domain2.end());
 	domain2.resize(distance(domain2.begin(),it2));
 	try{
 		if(domain.size() != domain2.size()) throw string(" Warning: Subunits labels are not unique in the .pdb file. Hope this is ok. ");
 	}
 	catch(const string & s){
 		cout << s << endl;
 	}
 	domain.clear();
 	domain=domain2;
 	map<string,int> pop;
 	for(size_t i=0;i<domain.size();i++)
 		pop[domain[i]]=i;

 	for(size_t i=0;i<MySub.size();i++){
 		string ss=MySub[i].getC();
 		MySub[i].putN(pop[ss]);
 		if(bAuto) MySub[i].putB(MyGroup.getBeta(i));
 	}
 	int ocount=0;

 	if(bAuto) {
 		for(size_t i=0;i<data.Size();i++){
 			string ss;
 			if(MySub[ocount].getB()) ss=" 1.00";
 			else ss=" 0.00";
 			*fDomainPDB << data[i].first << ss << data[i].last << endl;
 			ocount++;
 			//else *fDomainPDB  << data[i] << endl;
 		}
 	}

 // Interdomain stuff ends

  	unsigned int typesize=cidx0.size();
 	cidx=vector<vector<int> >(typesize);
 	TypeNames=vector<string>(typesize);

 	map<int,vector<int> >::iterator it=cidx0.begin();
 	for(int n=0;it != cidx0.end();++it, n++){
 		cidx[n]=it->second;
 	}
 	for(size_t n=0;n<cidx0.size();n++)
 		TypeNames[n]=mytype.strFind(n);

 	cidx0.clear();
 	std::sort(atres.begin(),atres.end(),op_sort());
 	std::transform(atres.begin(),atres.end(),atres.begin(),op_compare());

 	for(unsigned int i=0;i<atres.size();i++){
 		resind.push_back(atres[i].n);
  	}
 	vector<MyResidue>::iterator ita=std::unique(atres.begin(),atres.end(),myunique);
 	atres.resize(ita-atres.begin());

 	for(unsigned int i=0;i<atres.size();i++){
 		resinfo.push_back(atres[i].name);
 	}

 	nr=atomname.size();
 	nres=atres.size();
 	atResType=vector<int>(nr);

 	cout << std::fixed << std::setw(6)<< "\n Topology extracted from pdb file. Number of atoms: "
 			<< nr << " Number of residues: " << nres << "\n\n";
 	cout << std::showpoint << std::fixed << std::left;
 	cout << std::setw(17) << "  Residue Name"<< std::setw(16) << "No. of Atoms" << std::setw(17)
 	<< "Residue Type" << std::setw(17) << "Residue Type Numb." << endl;;

 	size_t AtomTotal=0;
 	for(unsigned int n=0;n<TypeNames.size();n++){
 		cout << std::left << std::fixed<< "    " <<std::setw(10) <<
 				TypeNames[n] << "      " << std::setw(12) << cidx[n].size()
 				<< "   " << std::setw(19) << mytype.strFindResType(n) << std::setw(17) << mytype.strFindResTypeNum(n)  << endl;
 		restype_s.push_back(mytype.strFindResType(n));
 		restype_i.push_back(mytype.strFindResTypeNum(n));
 		AtomTotal+=cidx[n].size();
 	}
 	cout <<"\n" << std::left << std::fixed<< "    " <<std::setw(10) << " Total  = " << "      "
 			<< std::setw(12) << AtomTotal << "\n" <<endl;
 	ResType=restype_i;
 	sort(ResType.begin(),ResType.end());
 	vector<int>::iterator itu=unique(ResType.begin(),ResType.end());
 	ResType.resize(std::distance(ResType.begin(),itu));

 	maxtype=TypeNames.size();

 	for(;maxtype>0;){
 		if(cidx[maxtype-1].size()) break;
 		maxtype--;
 	}
 	int o{0};
 	for(auto it: cidx){
 		for(auto jt: it){
 			int i=jt;
 			atResType[i]=o;
 		}
 		o++;
 	}
// 	FindGroups(xc,mytype.gWaterId());
}
} /* namespace Topol_NS */
