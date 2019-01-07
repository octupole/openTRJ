/*
 * FStream.h
 *
 *  Created on: May 21, 2012
 *      Author: marchi
 */

#ifndef FSTREAM_H_
#define FSTREAM_H_
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
using std::string;
using std::vector;
using std::ios;
using std::ifstream;
using std::cout;
using std::endl;
#define _FILE_OFFSET_BITS 64
const int HDIM=4;
const int TDIM=80;
const int FORTRANBYTES=4;

class Fstream {
protected:
	int dStep;
	int offStep;
	int nframe;
	int framenumber{0};
	int natoms{0};
	bool timescale{true};
	virtual void CompFrameNumber()=0;

public:
	Fstream():dStep(1),offStep(0),nframe(0){};
	Fstream(string & x): dStep(1),offStep(0),nframe(0){};
	virtual int gFrameStep()=0;
	virtual int gFrameNumber(){CompFrameNumber();return framenumber;}
	void setNatoms(int n){natoms=n;}
	bool CheckFrameNumber(int y){
		try{
			if(framenumber <= y) throw string("Number of frames is smaller than expected:");
		} catch(const string & s){
			cout << s << " it is " << framenumber-1 << " rather than " << y << endl;
			return false;
		}
		return true;
	}
	bool gTimescale(){return timescale;}
	int gStep(){return dStep;}
	int goffStep(){return offStep;}
	int gFrame(){return nframe;};
	virtual void Rewind()=0;
	virtual void seekg(off_t ,string)=0;
	virtual void nextFrame()=0;
	virtual ios::pos_type tellg()=0;
	virtual ~Fstream();
};

#endif /* FSTREAM_H_ */
