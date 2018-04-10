/*
 * ExecuteVoronoi.h
 *
 *  Created on: Nov 16, 2017
 *      Author: marchi
 */

#ifndef EXECUTEVORONOI_H_
#define EXECUTEVORONOI_H_
#include <cmath>
#include "FstreamC.h"
#include "FstreamF.h"
#include "MyUtilClass.h"
#include <iterator>

#include "../Execute/TrjRead.h"
#include "Finalize.h"
#include "Topol.h"
#include "FComms.h"
#include "IteratorAtoms.h"
#include "IteratorVoronoi.h"
#include "Finalize.h"
#include "Percolation.h"
#include "Atoms.h"
template<typename T>
class Atoms;

using namespace MATRIX;
using namespace DVECT;
using namespace Topol_NS;
namespace Voro {
class Voronoi;

template<typename T>
class ExecuteVoronoi {
	Topol * Top{nullptr}; ///< An instantiated Topol class pointer.
	Topol * TopPBC{nullptr};

	long int nstart{0},nend{-1},nskip{1}; ///< Where to start to end and how many steps to skip in between
	ios::streampos len;
	static size_t nnx,nny,nnz; ///< The three dimension of the grid
	bool Clustering{false}; ///< Do clustering or not
	 bool bTest{false};
	/// @cond TEST
	ofstream * foutp1{nullptr};
	int nacc{1};
	/// @endcond
	static Parallel::NewMPI * CurrMPI; ///< Point to current communicator if mpi invoked
	Parallel::FComms * Comms{nullptr};
	vector<size_t> MyRef;

	/// @cond TEST
	bool bDel{false};
	bool bHyd{false};
	bool binOutput{false};
	bool JSONOutput{false};
	double Rcut_in{-1.0},Rcut{-1.0};
	/// @endcond
	double MassSolute{-1.0}; ///< Molecular mass of the solute. Used to compute alpha.

	string fileout{""}; ///< filename of the ouput
	string fileout_bin{""}; ///< filename of the binary ouput

	ifstream * fpdb{nullptr}; ///< input stream pointer for pdb file
	Fstream * finx{nullptr}; ///< Interface for special input stream for trajectory files
	ofstream * foutx{nullptr}; ///< output stream
	ifstream * fin1x{nullptr}; ///<input stream pointer to saved voronoi binary file
	/// @cond TEST
	ifstream * fidb{nullptr};
	Voronoi * vor{nullptr};
	/** \brief
	 * private method to print the class
	 * @param y the output stream
	 */
	virtual void __Print(ostream & y);
	bool bOnce{false};

	/**
	 * Setup data members
	 * @param y the instantiated input class
	 */
	void __SetUp(trj::TrjRead & y);
	void __RunPDB(Atoms<T> *);
	void __RunTrajectory(Atoms<T> *);
	void __RunPost();
public:
	ExecuteVoronoi(trj::TrjRead & x);
	ExecuteVoronoi(trj::TrjRead & x, Topol & y);
	ExecuteVoronoi()=delete;
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
	/** \brief Initiate execution
	 *
	 * @param a MAtoms pointer
	 */
	void operator()(Atoms<T> *);


	/** \brief Given an ofstream print the class data
	 *
	 * @param fout ostream reference to ouput file
	 * @param y an instantiated Execbsaxs class
	 * @return
	 */
	friend ostream & operator<<(ostream & fout, ExecuteVoronoi & y){
		y.__Print(fout);
		return fout;
	}
	virtual ~ExecuteVoronoi();

};

} /* namespace Voronoi */

#endif /* EXECUTEVORONOI_H_ */
