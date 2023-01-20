/*
 * ExecbSaxs.h
 *
 *  Created on: Dec 16, 2015
 *      Author: marchi
 */

#ifndef SRC_EXECBSAXS_H_
#define SRC_EXECBSAXS_H_
#include <cmath>
#include "FstreamC.h"
#include "FstreamF.h"
#include "MyUtilClass.h"
#include "IteratorAtoms.h"
#include <iterator>
#include "NewMPI.h"
#include "Saxs.h"
#include "SaxsBSP.h"
#include "SaxsBSPStatic.h"
#include "SaxsBSPfixed.h"
#include "SaxsDebye.h"
#include "SaxsDir.h"
#include "RhoSaxs.h"
#include "RhoSaxsLI.h"
#include "RhoSaxsBSP.h"
#include "LagrangeInterpolation.h"
#include "PickSelection.h"
#include "TrjRead.h"
#include "fftw++.h"
#include "NewMPI.h"
#include "Finalize.h"
#include "DensMode.h"

using namespace MATRIX;
using namespace DVECT;
/*! \brief Base class and interface for
 *         classes running saxs calculations.
 *
 *  Detailed description starts here.
 */
class ExecbSaxs {
protected:
	using MAtoms=Atoms<double>;
	using Matrix=MATRIX::MMatrix<double>;
	using Dvect=DVECT::DDvect<double>;
	int nstart{0},nend{-1},nskip{1}; ///< Where to start to end and how many steps to skip in between
	std::streampos len;
	static size_t nnx,nny,nnz; ///< The three dimension of the grid
	bool bBSP{false};  ///< calculation with or without cardinal B-spline
	bool bFluct{true};  ///< Do the calculation including fluctuations
	bool bFixed{false};  ///< Do the calculation including fluctuations
	bool Clustering{false}; ///< Do clustering or not
	DensMode myDens;
	int myDensAvg{0};

	/// @cond TEST
	ofstream * foutp1{nullptr};
	int nacc{1};
	static double ddx1;
	static double ddx2;
	static double cutX;
	static double cutYZ;
	/// @endcond
	static int MyOrder; ///< Order of the saxs method
	static Parallel::NewMPI * CurrMPI; ///< Point to current communicator if mpi invoked

	/// @cond TEST
	ifstream * fin1x{nullptr};
	ifstream * fin2x{nullptr};
	/// @endcond
	vector<size_t> MyRef;

	/// @cond TEST
	bool bDel{false};
	bool bHyd{false};
	bool bIsrd{false};
	bool bSaxs{false};
	bool bEdens{false};
	bool bDensR{true};
	bool bDensQ{false};
	Enums::Compute WhatToDo{Enums::SAXS};

	/// @endcond
	bool bContrast{false}; ///< contrast is present or not
	bool bSaxsBSP{false};	///< Do SAXS with cardinal B-spline
	bool bDebye{false}; ///< Debey full calculation
	bool bDirect{false}; ///< Debey histogram formula
	bool bnoSplineOut{false}; ///< do not use spline when computing the profile
	double Mycut{4.0}; ///< Histogram cutoff in A-1
	double Myd{0.1}; ///< Histogram bin
	/// @cond TEST
	double Rcut{0};
	double Rcut_in{0};
	/// @endcond
	double MassSolute{-1.0}; ///< Molecular mass of the solute. Used to compute alpha.

	string fileout{""}; ///< filename of the ouput
	string fileout_bin{""}; ///< filename of the binary ouput

	ifstream * fpdb{nullptr}; ///< input stream pointer for pdb file
	Fstream * finx{nullptr}; ///< Interface for special input stream for trajectory files
	ofstream * foutx{nullptr}; ///< output stream
	ofstream * fout_binx{nullptr};///< binary output stream

	/// @cond TEST
	ifstream * fin_contrast{nullptr};
	ofstream * fout_saxsx{nullptr};
	/// @endcond
	Saxs * MySaxs{nullptr}; ///< Pointer to the Saxs class to be instantiated

	/** \brief
	 * private method to print the class
	 * @param y the output stream
	 */
	virtual void __Print(ostream & y);
	/**
	 * Setup data members
	 * @param y the instantiated input class
	 */
	void __SetUp(trj::TrjRead & y);
public:
/**
 * Create instance from input class
 * @param y input class TrjRead
 */
	ExecbSaxs(trj::TrjRead & y);
	ExecbSaxs()=delete;
	/** \brief Initiate communication class
	 *
	 */
	static void InitComm(){
		CurrMPI=new Parallel::NewMPI();
	}
	/** \brief Get current communicators for mpi runs
	 *
	 * @return  communication pointer
	 */
	static Parallel::NewMPI * CurrCom(){return CurrMPI;}
	virtual void bPrint(ofstream & )=0;
	/** \brief Initiate execution
	 *
	 * @param a MAtoms pointer
	 */
	virtual void operator()(MAtoms * a){};


	virtual ~ExecbSaxs(){};
	/** \brief Given an ofstream print the class data
	 *
	 * @param fout ostream reference to ouput file
	 * @param y an instantiated Execbsaxs class
	 * @return
	 */
	friend ostream & operator<<(ostream & fout, ExecbSaxs & y){
		y.__Print(fout);
		return fout;
	}
};


#endif /* SRC_EXECBSAXS_H_ */
