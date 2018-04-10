/*
 * FabComms.h
 *
 *  Created on: May 14, 2012
 *      Author: marchi
 */

#ifndef FCOMMS_H_
#define FCOMMS_H_
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iterator>

#include "../libtraj/NewMPI.h"
//#include "Communicator.hpp"

using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::fstream;

namespace Parallel {
template <typename T>
inline char * as_byte(T & y){
	return reinterpret_cast<char *> (&y);
}
class FComms {
	ofstream * fout{nullptr};
	string * filename{nullptr};
	string fileout;
	int nstart{0};
	int nend{-1};
	NewMPI * comm{nullptr};
public:
	FComms(NewMPI * CurrMPI){comm=CurrMPI;};
	FComms(NewMPI *,ofstream &, string &, int,int,int,int);
	FComms(NewMPI *,int,int,int,int);
	void Fopen(ofstream & foutx,string & fileoutx){
		fileout=fileoutx;
		foutx.open(fileoutx.c_str(),ios::trunc|ios::out|ios::in|ios::binary);
		fout=&foutx;
	}
	int getStart(){return nstart;}

	int getEnd(){return nend;}
	const string getFilename(int n){return filename[n];}
	ofstream & getStream(){return *fout;}
	void removeFiles(){
		if(comm->Get_Size() == 1) return;
		if(!comm->Get_Rank()) {
			for(size_t n=0;n<comm->Get_Size();n++){
				std::cout << "Deleting scratch file : " << filename[n] << std::endl;
				try{if(remove(filename[n].c_str()) != 0) throw string("---->> Temporary files cannot be removed! <<-----");}
				catch(const string & s){std::cout << s << std::endl;}
			}
		}
		comm->Barrier();
	}
	void closeStream(){fout->close();}
	void appendStreams();
	virtual ~FComms();
};

} /* namespace Parallel */
#endif /* FABCOMMS_H_ */
