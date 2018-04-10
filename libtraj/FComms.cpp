/*
 * FabComms.cpp
 *
 *  Created on: May 14, 2012
 *      Author: marchi
 */

#include "../libtraj/FComms.h"

namespace Parallel {

FComms::FComms(NewMPI * currMPI,ofstream & foutx,string & fileoutx,int nstartx,int nendx,int ntot,int nskip){
	comm=currMPI;
	int nend_x=(nendx == -1)?ntot-1:nendx;
	fileout=fileoutx;
	if(comm->Get_Size() > 1){
		int nrank=comm->Get_Rank();
		int nsize=comm->Get_Size();
		try{if(nsize > (nend_x-nstartx+1)) throw string("Number of processors is larger than the number of steps."
				" You are using more processors than you need!");}
		catch(const string & s){std::cout << s << std::endl;exit(-1);}
		try{if((nend_x-nstartx+1)%nsize) throw string("Warning number of steps is not divisible by the number of processors");}
		catch(const string & s){std::cout << s << std::endl;}

		int Myntot=(nend_x-nstartx+1)/nsize;
		filename=new string [nsize];
		nstart=nstartx+nrank*Myntot;
		nend=nstart+Myntot-1;
		vector<int> flen(nsize,0);
		if(!comm->Get_Rank()){
			for(int n=0;n<nsize;n++){
				char * buffer0= strdup("./.tmpfileXXXXXX");
				try{if(mkstemp(buffer0) == -1) throw "temporary file cannot be created";}
				catch(const char * s){std::cout << s << std::endl;exit(-1);}
				filename[n].assign(buffer0);
			}
			for(int n=0;n<nsize;n++) flen[n]=filename[n].size();
		}
		comm->Barrier();
		comm->Broadcast(&flen[0],nsize);
		int totflen=0;
		for(int n=0;n<nsize;n++) totflen+=flen[n];

		string names(totflen,' ');
		if(!comm->Get_Rank()){
			names.clear();
			for(int n=0;n<nsize;n++) names+=filename[n];
		}

		comm->Barrier();
		comm->Broadcast(&names[0],totflen);
		int pos=0;
		for(int n=0;n<nsize;n++) {
			filename[n].assign(names.substr(pos,flen[n]));
			pos+=flen[n];
		}
		if(foutx.is_open()) foutx.close();
		foutx.open(filename[nrank].c_str(),ios::trunc|ios::out|ios::in|ios::binary);
		fout=&foutx;
	} else{
		nstart=nstartx;
		nend=nendx;
		foutx.open(fileoutx.c_str(),ios::trunc|ios::out|ios::binary);
		fout=&foutx;
	}

}
FComms::FComms(NewMPI * currMPI,int nstartx,int nendx,int ntot,int nskip){
	comm=currMPI;

	int nend_x=(nendx == -1)?ntot-1:nendx;
	if(comm->Get_Size() > 1){
		int nrank=comm->Get_Rank();
		int nsize=comm->Get_Size();
		try{if(nsize > (nend_x-nstartx+1)) throw string("Number of processors is larger than the number of steps."
				" You are using more processors than you need!");}
		catch(const string & s){std::cout << s << std::endl;exit(-1);}
		try{if((nend_x-nstartx+1)%nsize) throw string("Warning number of steps is not divisible by the number of processors");}
		catch(const string & s){std::cout << s << std::endl;}

		int Myntot=(nend_x-nstartx+1)/nsize;
		nstart=nstartx+nrank*Myntot;
		nend=nstart+Myntot-1;
	} else{
		nstart=nstartx;
		nend=nendx;
	}
}
void FComms::appendStreams(){
	int nsize=comm->Get_Size();
	if(nsize < 2) return;
	int nrank=comm->Get_Rank();
	fout->close();
	comm->Barrier();
	if(!nrank) {
		fstream myfout(fileout.c_str(),fstream::out|fstream::trunc|fstream::binary);
		try{if(!myfout.is_open()) throw string("What??");}
		catch(const string & s){std::cout << s <<"\n";exit(-1);}

		for(int n=0;n<nsize;n++){
			fstream myfin(filename[n].c_str(),fstream::in|fstream::binary);
			myfin.seekg(0,ios::beg);
			myfout << myfin.rdbuf();
			myfin.close();
		}
		myfout.flush();
		myfout.close();
	}
	comm->Barrier();
}
FComms::~FComms() {}

} /* namespace Parallel */
