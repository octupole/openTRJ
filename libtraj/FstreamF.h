/*
 * FstreamF.h
 *
 *  Created on: May 21, 2012
 *      Author: marchi
 */

#ifndef FSTREAMF_H_
#define FSTREAMF_H_
#include "../libtraj/Fstream.h"

class FstreamF: public Fstream {
	ifstream fin;
public:
	FstreamF(string & file){fin.open(file.c_str(),ios::in|ios::binary);};
	ifstream & getfin(){return fin;}
	virtual bool CheckFrameNumber(int y){return true;};
	virtual void nextFrame();
	virtual void seekg(off_t,string);
	virtual ios::streampos tellg();
	virtual void Rewind(){
		nframe=0;
	}

	virtual ~FstreamF();
};

#endif /* FSTREAMF_H_ */
