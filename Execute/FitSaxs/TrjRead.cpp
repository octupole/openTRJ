/*
 * TrjRead.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: marchi
 */

#include "TrjRead.h"
namespace trj {



TrjRead::TrjRead(int nv,char ** v): trjInput::trjInput(nv,v) {
	// TODO Auto-generated constructor stub
	string comm0(v[0]);
	size_t mypos=comm0.find_last_of("/")+1;
	size_t length=comm0.size()-mypos;
	string command=comm0.substr(mypos,length);
	string errmsg;
	string Usage="Usage:\t"+ command + "\n";
	vector<string> use=this->getUsage();
	vector<string> SelRes;
	string Reference;
	for(unsigned int n=0;n<use.size();n++)
		Usage+=use[n];
	Usage+="\n\t Default values in square brackets []\n";
	try{
		if(nv == 1) throw Usage;
		if(int m=this->bTest().size()) {
			errmsg=" Command(s) not found: ";
			for(unsigned int n=0;n<m;n++)
				errmsg+=this->bTest()[n]+"  ";
			errmsg+="\n"+Usage;
			throw errmsg;
			}
		}
	catch(const string & s){
		cout << s << endl;
		exit(0);
	}

}
void TrjRead::Input(){
	try{
		if(!inmap["-in"].empty()) {
			if(inmap["-in"].size() < 3) throw string("\n two filenames expected for " + inmap["-in"][0] + " option\n ");
			if(inmap["-in"].size() > 3) throw string("\n More than two entries for " + inmap["-in"][0] + " option \n");
			filein1=inmap["-in"][1];
			ftest.open(filein1.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open function1 " + filein1 + "!!\n");
			ftest.close();
			fin1=new ifstream(filein1);
			filein2=inmap["-in"][2];
			ftest.open(filein2.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open function2 " + filein2 + "!!\n");
			ftest.close();
			fin2=new ifstream(filein2);
			cout << filein1 << " " << filein2 <<endl;
		}
		if(!inmap["-notr"].empty()) {
			if(inmap["-notr"].size() > 1) throw string(inmap["-notr"][0] + " option need no argument!\n");
			noTransl=true;
		}
		if(!inmap["-o"].empty()) {
			if(inmap["-o"].size() < 2) throw string("\n filename expected for " + inmap["-o"][0] + " option \n");
			if(inmap["-o"].size() > 2) throw string("\n More than one entry for " + inmap["-o"][0] + " option \n");
			fileout=inmap["-o"][1];
		}
	}catch(const string & s){
		cout << s << endl;
		exit(0);
	}

	foutx=new ofstream(fileout.c_str(),ios::out);


	gFin1=fin1;
	gFin2=fin2;
	gFoutx=foutx;
	try{
		if(!(this->gFin1() && this->gFin2())) throw string(" No input file(s) given program stops! ");
	} catch(const string & s){
		cout << s <<endl;
		exit(0);
	}
}
TrjRead::~TrjRead() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */
