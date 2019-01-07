/*
 * FstreamF.cpp
 *
 *  Created on: May 21, 2012
 *      Author: marchi
 */

#include "../libtraj/FstreamF.h"

void FstreamF::seekg(off_t n,string x){
	try{
		if(x.find("beg") != string::npos)
			fin.seekg(static_cast<ios::streamoff> (n),ios::beg);
		else if(x.find("end") != string::npos)
			fin.seekg(static_cast<ios::streamoff> (n),ios::end);
		else
			throw string("Seek direction unknown ");
	}
	catch(const string & s){
		cout << s << endl;
		exit(-1);
	}
}
ios::pos_type FstreamF::tellg(){
	return fin.tellg();
}
FstreamF::FstreamF(string & file){
	fin.open(file.c_str(),ios::in|ios::binary);
	CompFrameNumber();
	timescale=false;
}

void FstreamF::Rewind(){
	fin.seekg(headerPos,ios::beg);
}
void FstreamF::CompFrameNumber(){
	char hdr[HDIM];
	int nfr{0},istart{0},nsavc{0},natoma{0},tmp0[5],tmp1[9];
	double delta{0};
	int ntitle;
	char tmp[TDIM];
	vector<string> title;

	fin.seekg(FORTRANBYTES,ios::beg);
	fin.read(reinterpret_cast<char *> (&hdr),sizeof(hdr));
	fin.read(reinterpret_cast<char *> (&nfr),sizeof(nfr));
	fin.read(reinterpret_cast<char *> (&istart),sizeof(istart));

	fin.read(reinterpret_cast<char *> (&nsavc),sizeof(nsavc));
	fin.read(reinterpret_cast<char *> (&tmp0),sizeof(tmp0));
	fin.read(reinterpret_cast<char *> (&natoma),sizeof(natoma));
	fin.read(reinterpret_cast<char *> (&delta),sizeof(delta));
	fin.read(reinterpret_cast<char *> (&tmp1),sizeof(tmp1));
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.read(reinterpret_cast<char *> (&ntitle),sizeof(ntitle));

	for(int i{0};i< ntitle;i++){
		fin.read(reinterpret_cast<char *> (&tmp),sizeof(tmp));
	}
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.read(reinterpret_cast<char *> (&natoms),sizeof(natoms));
	fin.seekg(FORTRANBYTES,ios::cur);

	headerPos= fin.tellg();
	fin.seekg(0,ios::end);
	off_t totTrj=fin.tellg()-headerPos;

	fin.seekg(0,ios::beg);
	off_t totRecord=FORTRANBYTES*2+FORTRANBYTES*6+sizeof(double[6])+sizeof(float)*3*natoms;
	framenumber= totTrj/totRecord;
}

void FstreamF::nextFrame(){
	nframe+=dStep;
}

FstreamF::~FstreamF() {
	// TODO Auto-generated destructor stub
}

