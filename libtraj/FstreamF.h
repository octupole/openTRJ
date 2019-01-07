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
	off_t headerPos{0};
public:
	FstreamF(string & file);
	ifstream & getfin(){return fin;}
	int gFrameStep(){return framenumber;}
	void CompFrameNumber();
	void nextFrame();
	void seekg(off_t,string);
	ios::streampos tellg();
	void Rewind();

	virtual ~FstreamF();
};

#endif /* FSTREAMF_H_ */
