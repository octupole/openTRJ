/*
 * FstreamC.cpp
 *
 *  Created on: May 21, 2012
 *      Author: marchi
 */

#include "../libtraj/FstreamC.h"
FstreamC::FstreamC(string file,string s){
	if(s.substr(0,1)=="r"){
		fin = xdrfile_open(file.c_str(),"r");
		char * myFile=new char[file.size()];
		strcpy(myFile, file.c_str());
		read_xtc_natoms(myFile,&natoms);
		sStep();
		CompFrameNumber();
	}else{
		fin = xdrfile_open(file.c_str(),"w");
	}

}
void FstreamC::seekg(off_t n,string x){
	try{
		FILE * fp=xdrfile_get_fp(fin);
		if(x.find("beg") != string::npos)
			rewind(fp);
		else if(x.find("end") != string::npos)
			fseeko(fp,n,SEEK_END);
		else if(x.find("cur") != string::npos)
			fseeko(fp,n,SEEK_CUR);
		else
			throw string("Seek direction unknown ");
	}
	catch(const string & s){
		cout << s << endl;
		exit(-1);
	}
}
int FstreamC::gFrameStep(){
	FILE * fp=xdrfile_get_fp(fin);
	XDR * myxdr=xdrfile_get_xdr(fin);
	ios::pos_type myCurr=tellg();
	rewind(fp);
	int f0=xtc_get_next_frame_number(fp,myxdr,natoms);
	fseeko(fp,0,myCurr);
	return (framenumber-f0)/dStep;
}


void FstreamC::CompFrameNumber(){
	FILE * fp=xdrfile_get_fp(fin);
	XDR * myxdr=xdrfile_get_xdr(fin);

	int bOK{0};
	int f0=xtc_get_next_frame_number(fp,myxdr,natoms);

	framenumber= xdr_xtc_get_last_frame_number(fp,myxdr,natoms,&bOK);

}

void FstreamC::nextFrame(){
	nframe+=dStep;
}

void FstreamC::sStep(){
	rvec   *x;
	matrix box;
	int    step0;
	float   prec,time0;

	FILE * fp=xdrfile_get_fp(fin);
	XDR * myxdr=xdrfile_get_xdr(fin);
	rewind(fp);
	int f1;
	try{
		f1=xtc_get_next_frame_number(fp,myxdr,natoms);
		x=new rvec[natoms];
		if(read_xtc(fin,natoms,&step0,&time0,box,x,&prec))
			throw string("Cannot read first frame!");
		rewind(fp);
		delete [] x;
	}
	catch(const string & s){
		cout << s << endl;
		exit(-1);
	}
	dStep=f1-step0;
	offStep=step0;
	nframe=step0;
};

ios::pos_type FstreamC::tellg(){
	FILE * fp=xdrfile_get_fp(fin);
	return ftello(fp);
}

FstreamC::~FstreamC() {
	// TODO Auto-generated destructor stub
}

