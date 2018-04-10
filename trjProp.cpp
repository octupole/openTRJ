//============================================================================
// Name        : trjProp.cpp
// Author      : Massimo Marchi
// Version     : 0.2alpha
// Copyright   : CeCILL copyright
//============================================================================
/*! \mainpage trjVoronoi
 *
 * \section intro_sec Overview
 *
 *
 * \section install_sec Installation
 *
 *	trjVornoi uses a standard autoconf script and makefiles created by
 *	automake, like most GNU programs. This means your normal installation
 *	actions will be limited to run a
 *	sh bootstrap (generates all autoconf file needed)
 *
 *	./configure
 *	make
 *	make install
 *
 *	The code is dependent on the following packages:
 *
 *	i)  A modified xdrfile library
 *	ii) The voro++ library.
 *
 *	Both are provided with this distribution, but voro++ is integrated
 *	with the main code, whereas xdrfile need to be compiled separately.
 *
 *	Before installing trjVoronoi, go to the xdrfile directory and issue:
 *
 *	sh bootstrap
 *	./configure
 *	make install
 *
 *	Following the installation of the xdr package, return to the installation
 *	directory and issue a
 *	sh bootstrap
 *	./configure
 *	make
 *	make install
 *
 *	Noticeably ./configure has amongst others the following options:
 *
 *	--enable-intel		compile trjVoronoi without user-define literals to
 *					    compile on an intel compiler
 *	--enable-mpi		compile trjVoronoi with MPI support
 *
 *	trjVoronoi uses several features available with the C++11 standard, thus
 *	it will not compile with GNU g++ compilers earlier than 4.8 and with
 *	Intel icpc compilers versions earlier than 15.0.2. It has been
 *	installed on Centos and Debian Linux machines and on Mac OS X.
 */
#include <iostream>
#include <vector>
#include <string>

#include <iterator>
#include "TopolPDB.h"
#include "Topol.h"
#include "Timer.h"
#include "myEnums.hpp"
#include "Timer.h"
#include "NewMPI.h"
#include "Atoms.h"
#include "Execute/ClearUsage.h"
#include "Execute/ExecuteProp.h"
#include "Execute/TrjRead.h"
#include "PickSelection.h"
#include "Finalize.h"

using namespace Topol_NS;

using namespace std;

Atoms<double> * atm{nullptr};
Properties::ExecuteProp<double> * MyRun{nullptr};


TopolPDB topPDB;

#include <iostream>
using namespace std;
int main(int argc, char ** argv){
	/*
	**
	 * set communicator
	 */
	Properties::ExecuteProp<double>::InitComm();// Start parallel tasks
	Finale::Finalize::CurrMPI=Properties::ExecuteProp<double>::CurrCom();
	trj::TrjRead::SetComm(Properties::ExecuteProp<double>::CurrCom());
	/**
	 * read input
	 */
	ClearUsage clr({6,7,9,10,13,14,16,18,19,20,21,22,23,24});
	trj::TrjRead MyIn(argc,argv,clr);
	MyIn.Input();
	timer::Timer myTime;
	Topol_NS::Topol MyTop;

	if(MyIn.gFin1()){
		MyRun=new Properties::ExecuteProp<double>(MyIn);
	} else {
		vector<string> data;
		if(MyIn.gFpdb()){
			// read pdb file to construct topology
			for(string str;getline(*MyIn.gFpdb(),str);){
				data.push_back(str);
			}
		}
		topPDB.sDetPolsegment(MyIn.gsDetResidue(),MyIn.gsPolarAtoms());
		topPDB(data);
		(MyTop)(topPDB,MyIn.bbIsrd());
		typedef map<string,vector<vector<int> > > mappa;
		mappa & Def=MyTop.getDef();
		MyIn.gReference()=PickSelection(MyIn.gReference()).Select<Enums::Reference>(Def,MyTop); // Pick the Reference residue

		int natoms=MyTop.Size();
		atm=new Atoms<double>(natoms);

		MyTop.InitSelection(MyIn.gReference(),Enums::Reference);
		MyRun=new Properties::ExecuteProp<double>(MyIn,MyTop);

	}

	(*MyRun)(atm);
	double Steps=(double) (MyIn.gnend()-MyIn.gnstart())/(double) MyIn.gnskip();
	double nsteps=Steps<=0?1.0:Steps;
	cout << "\n----> Elapsed Time per step:  " << myTime/nsteps << " sec " <<endl;
	return 0;
}
