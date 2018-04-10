/*
 * ExecbSaxsPost.h
 *
 *  Created on: Mar 6, 2016
 *      Author: marchi
 */

#ifndef SRC_EXECBSAXSPOST_H_
#define SRC_EXECBSAXSPOST_H_

#include "ExecbSaxs.h"
/** \brief Generate instances of the saxs classes and and use saxs methods to obtain saxs profile
 *
 */
class ExecbSaxsPost: public ExecbSaxs {
	/** \brief Copy some of the TrjClass data to the class members.
	 *         All members are declared in the base class
	 *
	 * @param x Reference to the input class
	 */
	void __SetUp(trj::TrjRead & x);
	/** \brief Compute differences between Intensities, then write the histogram of the profile
	 *
	 */
	void __Differences();
	void __GofR();
public:
	ExecbSaxsPost()=delete;
	/** \brief The only constructor. Must be defined after the input class is instantiated
	 *
	 * @param x
	 */
	ExecbSaxsPost(trj::TrjRead & x);
	/** \brief Function call operator
	 *
	 * @param x a pointer to an MAtoms class. x is not used and it might be a nullptr.
	 *          Kept for compatibility with the base class.
	 */
	virtual void operator()(MAtoms * x);
	/** \brief binary print which does nothing in this class
	 *
	 * @param fout ofstream
	 */
	virtual void bPrint(ofstream & fout){};

	virtual ~ExecbSaxsPost();
};

#endif /* SRC_EXECBSAXSPOST_H_ */
