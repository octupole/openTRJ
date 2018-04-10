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

#include <Fstream.h>
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <xdrfile_seek.h>


class FstreamC: public Fstream {
	XDRFILE * fin{nullptr};
	void sStep();
	int framenumber{0};
	void CompFrameNumber();
	int natoms{0};
public:
	FstreamC(string ,string);
	virtual int gFrameNumber(){CompFrameNumber();return framenumber;}
	virtual int gFrameStep();

	XDRFILE * getfin(){return fin;}
	int gNatoms(){return natoms;}
	int gStep(){return dStep;}
	virtual bool CheckFrameNumber(int y){
		try{
			if(framenumber <= y) throw string("Number of frames is smaller than expected:");
		} catch(const string & s){
			cout << s << " it is " << framenumber-1 << " rather than " << y << endl;
			return false;
		}
		return true;
	}
	virtual void Rewind(){
		FILE * fp=xdrfile_get_fp(fin);
		rewind(fp);

	};

	virtual void nextFrame();
	virtual void seekg(off_t,string);

	virtual ios::streampos tellg();
	virtual ~FstreamC();
};

#endif /* FSTREAMC_H_ */
