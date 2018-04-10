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
using std::string;
using std::ios;
using std::ifstream;
using std::cout;
using std::endl;
#define _FILE_OFFSET_BITS 64

class Fstream {
protected:
	int dStep;
	int offStep;
	int nframe;
public:
	Fstream():dStep(1),offStep(0),nframe(0){};
	Fstream(string & x): dStep(1),offStep(0),nframe(0){};
	virtual int gFrameNumber(){return 0;}
	virtual int gFrameStep(){return 0;}

	int gStep(){return dStep;}
	int goffStep(){return offStep;}
	int gFrame(){return nframe;};
	virtual void Rewind()=0;
	virtual void seekg(off_t ,string)=0;
	virtual bool CheckFrameNumber(int)=0;
	virtual void nextFrame()=0;
	virtual ios::pos_type tellg()=0;
	virtual ~Fstream();
};

#endif /* FSTREAM_H_ */
