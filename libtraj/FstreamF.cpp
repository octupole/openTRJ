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

void FstreamF::nextFrame(){
	nframe+=dStep;
}

FstreamF::~FstreamF() {
	// TODO Auto-generated destructor stub
}

