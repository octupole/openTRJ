/*
 * FstreamC.h
 *
 *  Created on: May 21, 2012
 *      Author: marchi
 */

#ifndef FSTREAMC_H_
#define FSTREAMC_H_
#include <cstdlib>
#include <cstring>

#include "Fstream.h"
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <xdrfile_seek.h>


class FstreamC: public Fstream {
	XDRFILE * fin{nullptr};
	void sStep();
	int framenumber{0};
	void CompFrameNumber();
public:
	FstreamC(string ,string);
	virtual int gFrameNumber(){CompFrameNumber();return framenumber;}
	virtual int gFrameStep();

	XDRFILE * getfin(){return fin;}
	int gNatoms(){return natoms;}
	int gStep(){return dStep;}
	virtual void Rewind(){
		FILE * fp=xdrfile_get_fp(fin);
		rewind(fp);

	};

	virtual void nextFrame();
	virtual void seekg(off_t,string);

        virtual std::streampos tellg();
	virtual ~FstreamC();
};

#endif /* FSTREAMC_H_ */
