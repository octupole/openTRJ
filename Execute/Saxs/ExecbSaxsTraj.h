/*
 * ExecbSaxs.h
 *
 *  Created on: Dec 16, 2015
 *      Author: marchi
 */

#ifndef SRC_EXECBSAXSTRAJ_H_
#define SRC_EXECBSAXSTRAJ_H_
#include "ExecbSaxs.h"
#include "Contacts.h"
#include "myEnums.hpp"
#include "SaxsFilter.h"


using namespace MATRIX;
using namespace DVECT;


/** \brief Generate instances of the saxs classes and applies saxs methods to the trajectory
 *
 */
using MAtoms=Atoms<double>;
using Matrix=MATRIX::MMatrix<double>;
using Dvect=DVECT::DDvect<double>;

class ExecbSaxsTraj: public ExecbSaxs {
	using PercolationD=Percolation<double>;
	using MetricD=Metric<double>;
	using ContactsD=Contacts<double>;
protected:
	Topol * Top{nullptr}; ///< An instantiated Topol class pointer.
	RhoSaxs * Rho_ex{nullptr}; ///< A pointer for the SAXS density. The class instantiate Rho_ex to
							   ///< the wanted derived class of Rhosaxs
							   ///
	Topol * TopPBC{nullptr};
	SaxsFilter Filter;
	Enums::Padding exPadding{Enums::myDensity};
	double SuperCell0; ///< The supercell scaling factor
	double SuperCell; ///< The largest sides
	double Rcut_in{-1.0},Rcut{-1.0};
	bool bOnce{false};

	/** \brief Private method to iterate over the step of the trajectory
	 *
	 * @param x Instantiated MAtoms class
	 */
	void __RunTrajectory(MAtoms * x);
	/** \brief Private method to compute the SAXS of the coordinates of the pdb file
	 *
	 * @param x Instantiated MAtoms class
	 */

	void __RunPDB(MAtoms * x);
	void __RunPDBtest(MAtoms * x);

	/** \brief Private method setup the class member data based on input
	 *
	 * @param x Instantiated MAtoms class
	 */
	void __SetUp(trj::TrjRead & x);
	/** Set the supercell
	 *
	 * @param x Pointer to an instantiated MAtoms class
	 */
	void __SuperCell(MAtoms * x);
	/** \brief Allocate Rho_ex
	 *
	 * @param x An instantiated MAtoms class
	 */
	void __AllocateRho(MAtoms * x);

	/** \brief Make two dummy sites to be used to construct a the pbc saxs filter
	 *
	 */
	void __MakeDummyAtoms();

	void __Compute(const MAtoms *);

	MAtoms * Patm{nullptr};

public:
	/** \brief Constructor operator for the class. Needs two arguments not one like the base class.
	 *
	 * @param x Input class
	 * @param y Topology class
	 */
	ExecbSaxsTraj(trj::TrjRead & x, Topol & y);
	ExecbSaxsTraj()=delete;

	virtual ~ExecbSaxsTraj(){};
	/** \brief the call operator. This is a public method used to run a SAXS calculation
	 *
	 * @param x MAtoms pointer
	 */
	virtual void operator()(MAtoms * x);
	/** \brief Print the binary file
	 *
	 * @param x Print on ofstream x the class member data
	 */
	virtual void bPrint(ofstream & x);
	friend ofstream & operator<<(ofstream &, ExecbSaxsTraj &);
};


#endif /* SRC_EXECBSAXS_H_ */
