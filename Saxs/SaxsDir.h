/*
 * SaxsDir.h
 *
 *  Created on: Mar 19, 2016
 *      Author: marchi
 */

#ifndef TRJLIB_SAXSDIR_H_
#define TRJLIB_SAXSDIR_H_

#include "Saxs.h"
/** \brief Extends saxs to histogramatic Debey calculation through the ComputeSF method
 *
 */
class SaxsDir: public Saxs {
public:
	SaxsDir();
	SaxsDir(double dq0, double qcut0): Saxs(dq0,qcut0){};
	SaxsDir(int MyOrder,double dq0, double qcut0): Saxs(MyOrder,dq0,qcut0){};


	virtual void ComputeSAXS(RhoSaxs *,const MAtoms *);
	virtual ~SaxsDir();
};

#endif /* TRJLIB_SAXSDIR_H_ */
