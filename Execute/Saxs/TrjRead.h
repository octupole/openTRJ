/*
 * TrjRead.h
 *
 *  Created on: Dec 23, 2015
 *      Author: marchi
 */

#ifndef SRC_TRJREAD_H_
#define SRC_TRJREAD_H_

#include "trjInput.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iterator>
#include "FstreamC.h"
#include "FstreamF.h"
#include "CenterMass.h"
#include "Parameters.h"
#include "SaxsPadding.h"
#include "Finalize.h"
#include "BackupFile.h"
#include "NewMPI.h"
#include "Finalize.h"
#include "myEnums.hpp"
#include "DensMode.h"

#ifndef __INTEL
#include "Units.hh"
using namespace Units;
#endif

using namespace std;

namespace trj {
template <typename T>
class Streams{
	T * myStream{nullptr};
public:
	Streams(){};
	Streams(T * y):myStream{&y}{};
	 Streams & operator=(T * y ){
		myStream=y;
		return *this;
	}
	T * operator()(){
		return myStream;
	}

};

template <typename T>
class Values{
	T * const myValue;
public:
	Values(T & y): myValue{&y}{};
	T & operator()(){return *myValue;}
};

/** \brief The class extends trjInput by providing methods to read the on line options
 *         and provide the stored data to calling classes.
 *
 */
class TrjRead: public trjInput {
	using Dvect=DVECT::DDvect<double>;
	static Parallel::NewMPI * CurrMPI;
	int nstart{0};
	int nend{-1};
	int nskip{1};
	int nacc{1};
	// input read
	bool inputfile{false};
	bool bnoSplineOut{false};
	bool bFluct{true};
	bool bFixed{false};
	string filein,fileout="Micelles.out",filechg="Charges.dat",filepdb;
	string fileoutp1,fileoutp2,fileoutp3;
	string fileout_bin{"Micelles.bin"};


	ifstream fpdb,fchg;
	ifstream ftest,fdefdomain;
	ofstream fdomain;
	Fstream * finx{nullptr};
	Fstream * fin2x{nullptr};
	Fstream * fout_xtcx{nullptr};
	ofstream * fout_pdbx{nullptr};
	ofstream * fout_dip{nullptr};
	ofstream * fout_cmx{nullptr};
	ofstream * fout_elx{nullptr};
	ofstream * foutx{nullptr};
	ofstream * fout_binx{nullptr};
	ofstream * foutp1{nullptr};
	ofstream * foutp2{nullptr};
	ofstream * foutp3{nullptr};
	ofstream * fclust{nullptr};
	ifstream * fin1{nullptr};
	ifstream * fin2{nullptr};
	ifstream * fin_rcmx{nullptr};
	ifstream * fin_contrast{nullptr};
	ifstream * fin_padding{nullptr};
	ofstream * fout_saxsx{nullptr};
	vector<int> WignerArgs{0,0,0};
	Dvect e0{0.0};
	double Rcut{3.5};
	double Rcut_in{-1.0};
	double MassSolute{-1.0};
	double PercoCutoff{0.0};

	matrix MyCO={{0,0,0},{0,0,0},{0,0,0}};
	vector<string> SelRes;
	vector<string> Reference;
	SaxsPadding fftPadding;
	bool bDel{true};
	bool bWriteFab{false};
	bool bTestVol{false};
	bool bHyd{true};
	bool bIsrd{false};
	bool bVoro{false};
	bool bContrast{false};
	bool bRho{false};
	bool bprof{false};
	bool bClust{false};
	bool bOnce{false};
	bool bdip{false};
	bool bedens{false};
	bool bLagrange{false};
	bool bBSP{false};
	bool bOutBin{false};
	bool bNonEq{false};
	bool bSaxs{true};
	bool bDebye{false};
	bool bDirect{false};
	bool bSaxsBSP{false};
	bool bSans{false};
	bool bPost{false};
	bool bElDens{false};
	DensMode ModeCompute;

	int MyOrder{1};
	int MyDensAvg{4};
	size_t BoxMultiply{1};
	double SuperCellSide{1.0};
	double Mycut{4.0};
	double Myd{0.05};
	double ewSigma{0.00707};
	double alpha{0.0};
	double ddx1=0.5*unit_nm; // in Angstroems.
	double ddx2=1.5*unit_nm; // in Angstroems.
	double cutX{0.0};
	double cutYZ{0.0};
	unsigned int nnx{128},nny{128},nnz{128};
	string sVect="X";
	double ConvFactor{1.0};
public:
	Enums::Padding WhichPadding{Enums::avgDensity};
	Enums::Compute WhatToDo{Enums::SAXS};
	static void SetComm(Parallel::NewMPI * y){
		CurrMPI=y;
	}
	Values<vector<int> > gWignerArgs{WignerArgs};
	Values<Dvect> gE0{e0};

	Values<int> gnstart{nstart};
	Values<int> gnend{nend};
	Values<int> gnskip{nskip};
	Values<int> gnacc{nacc};
	Values<double> gRcut{Rcut};
	Values<double> gMassSolute{MassSolute};
	Values<double> gPercoCutoff{PercoCutoff};
	Values<double> gRcut_in{Rcut_in};
	Values<double> gFactor{ConvFactor};
	Values<DensMode> gModeCompute{ModeCompute};
	Values<vector<string> > gSelRes{SelRes};
	Values<vector<string> >gReference{Reference};
	Values<SaxsPadding>gfftPadding{fftPadding};
	CenterMass_t WhichDiffusion{diffk};

	Values<string> gfilein{filein};
	Values<string> gfileout{fileout};
	Values<string> gfileout_bin{fileout_bin};
	Values<string> gfilechg{filechg};
	Values<string> gfilepdb{filepdb};
	Values<string> gfileoutp1{fileoutp1};
	Values<string> gfileoutp2{fileoutp2};
	Values<string> gfileoutp3{fileoutp3};
	Values<string> gsVect{sVect};

	Values<size_t> gBoxMultiply{BoxMultiply};
	Values<double> gSuperCellSide{SuperCellSide};
	Values<int> gMyOrder{MyOrder};
	Values<int> gMyDensAvg{MyDensAvg};
	Values<unsigned int> gnnx{nnx};
	Values<unsigned int> gnny{nny};
	Values<unsigned int> gnnz{nny};
	Values<double> gddx1{ddx1};
	Values<double> gddx2{ddx2};
	Values<double> gcutX{cutX};
	Values<double> gcutYZ{cutYZ};
	Values<double> gMyCut{Mycut};
	Values<double> gMyd{Myd};
	Values<double> gewSigma{ewSigma};
	Values<double> galpha{alpha};
	Values<bool> binputfile{inputfile};
	Values<bool> bbnoSplineOut{bnoSplineOut};
	Values<bool> bbFluct{bFluct};
	Values<bool> bbFixed{bFixed};
	Values<bool> bbDel{bDel};
	Values<bool> bbWriteFab{bWriteFab};
	Values<bool> bbTestVol{bTestVol};
	Values<bool> bbHyd{bHyd};
	Values<bool> bbIsrd{bIsrd};
	Values<bool> bbVoro{bVoro};
	Values<bool> bbContrast{bContrast};
	Values<bool> bbRho{bRho};
	Values<bool> bbprof{bprof};
	Values<bool> bbClust{bClust};
	Values<bool> bbOnce{bOnce};
	Values<bool> bbdip{bdip};
	Values<bool> bbedens{bedens};
	Values<bool> bbLagrange{bLagrange};
	Values<bool> bbBSP{bBSP};
	Values<bool> bbOutBin{bOutBin};

	Values<bool> bbNonEq{bNonEq};
	Values<bool> bbSaxs{bSaxs};
	Values<bool> bbSans{bSans};
	Values<bool> bbDebye{bDebye};
	Values<bool> bbDirect{bDirect};
	Values<bool> bbPost{bPost};
	Values<bool> bbSaxsBSP{bSaxsBSP};


	Streams<Fstream> gFinx;
	Streams<Fstream> gFin2x;
	Streams<Fstream> gFout_xtcx;
	Streams<ofstream> gFout_pdbx;
	Streams<ofstream> gFout_dip;
	Streams<ofstream> gFout_cmx;
	Streams<ofstream> gFout_elx;
	Streams<ofstream> gFoutx;
	Streams<ofstream> gFout_binx;
	Streams<ofstream> gFoutp1;
	Streams<ofstream> gFoutp2;
	Streams<ofstream> gFoutp3;
	Streams<ofstream> gFclust;
	Streams<ofstream> gFoutsaxs;
	Streams<ifstream> gFincontrast;
	Streams<ofstream> gFdomain;
	Streams<ifstream> gFidb;
	Streams<ifstream> gFin1;
	Streams<ifstream> gFin2;
	Streams<ifstream> gFin_rcmx;
	Streams<ifstream> gFin_contrast;
	Streams<ifstream> gFin_padding;
	Streams<ifstream> gFpdb;
	Streams<ifstream> gFchg;
	Streams<ifstream> gFdefdomain;
	Streams<ifstream> gFtest;

	TrjRead(int nv,char ** v);
	void Input();
	virtual ~TrjRead();
};

} /* namespace trj */

#endif /* SRC_TRJREAD_H_ */
