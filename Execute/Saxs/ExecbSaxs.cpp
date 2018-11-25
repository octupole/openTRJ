/*
 * ExecbSaxs.cpp
 *
 *  Created on: Dec 16, 2015
 *      Author: marchi
 */

#include "ExecbSaxs.h"


double ExecbSaxs::ddx1=0.5*unit_nm; // in Angstroems
double ExecbSaxs::ddx2=1.5*unit_nm; // in Angstroems.
double ExecbSaxs::cutX=5.0*unit_nm;
double ExecbSaxs::cutYZ=40.0*unit_nm;
int ExecbSaxs::MyOrder=1;
size_t ExecbSaxs::nnx=1;
size_t ExecbSaxs::nny=1;
size_t ExecbSaxs::nnz=1;
Parallel::NewMPI * ExecbSaxs::CurrMPI{nullptr};
enum Model {Consecutive, Chunks};
const Model MyModel{Chunks};

ExecbSaxs::ExecbSaxs(trj::TrjRead & MyIn){

	__SetUp(MyIn);

};

/** \brief setup all the class members from the input handling class
 *
 * @param MyIn Input handling class trj::TrjRead
 */
void ExecbSaxs::__SetUp(trj::TrjRead & MyIn){
	ddx1=MyIn.gddx1();
	cutYZ=MyIn.gcutYZ();

	bDel=MyIn.bbDel();
	bHyd=MyIn.bbHyd();
	bIsrd=MyIn.bbIsrd();
	bContrast=MyIn.bbContrast();
	Mycut=MyIn.gMyCut();
	Myd=MyIn.gMyd();
	bEdens=MyIn.bbedens();

	fileout_bin=MyIn.gfileout_bin();
	fout_binx=MyIn.gFout_binx();

	fpdb=MyIn.gFpdb();
	finx=MyIn.gFinx();
	foutx=MyIn.gFoutx();

	fin_contrast=MyIn.gFin_contrast();
	fout_saxsx=MyIn.gFoutsaxs();
	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
	bnoSplineOut=MyIn.bbnoSplineOut();
	MassSolute=MyIn.gMassSolute();
	WhatToDo=MyIn.WhatToDo;

// To compute density

	myDens=MyIn.gModeCompute();
	myDensAvg=MyIn.gMyDensAvg();
	fileout=MyIn.gfileout();
}
void ExecbSaxs::__Print(ostream & fout ){
	if(!CurrMPI->Get_Rank()){
		fout << *MySaxs;
		cout << "\nProgram completed: Output data written to " + fileout << "\n\n";
	}

}

