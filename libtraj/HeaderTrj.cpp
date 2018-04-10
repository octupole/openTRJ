/*
 * HeaderTrj.cpp
 *
 *  Created on: May 13, 2012
 *      Author: marchi
 */

#include "../libtraj/HeaderTrj.h"

HeaderTrj::HeaderTrj() {
	// TODO Auto-generated constructor stub

}
void HeaderTrj::ReadHeader(FstreamC * fin){
	XDRFILE * xd=fin->getfin();
	rvec   *x;
	matrix box;
	int    step;
	float   prec,time;

	natoms=fin->gNatoms();
	x=new rvec[natoms];
	try{

		if(read_xtc(xd,natoms,&step,&time,box,x,&prec))
			throw string("Cannot read first frame!");

	}
	catch(const string & s){
		cout << s << endl;
		exit(-1);
	}
	nfr=step;
	istart=step;
	fin->Rewind();
	delete [] x;
}

void HeaderTrj::ReadHeader(FstreamF * fin){
	ReadHeader(fin->getfin());
}

void HeaderTrj::ReadHeader(ifstream & fin){
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.read(reinterpret_cast<char *> (&hdr),sizeof(hdr));
	fin.read(reinterpret_cast<char *> (&nfr),sizeof(nfr));
	fin.read(reinterpret_cast<char *> (&istart),sizeof(istart));
	const int OFFSET=92;
	fin.seekg(OFFSET,ios::beg);
	fin.seekg(FORTRANBYTES,ios::cur);
	int ntitle;
	fin.read(reinterpret_cast<char *> (&ntitle),sizeof(ntitle));
	for(int i=0;i<ntitle;i++){
		char tmp[TDIM];
		fin.read(reinterpret_cast<char *> (&tmp),sizeof(tmp));
		title.push_back(string(tmp));
	}
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.read(reinterpret_cast<char *> (&natoms),sizeof(natoms));
	fin.seekg(FORTRANBYTES,ios::cur);
}


HeaderTrj::~HeaderTrj() {
	// TODO Auto-generated destructor stub
}

