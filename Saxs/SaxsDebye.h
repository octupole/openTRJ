/*
 * SaxsDebey.h
 *
 *  Created on: Mar 19, 2016
 *      Author: marchi
 */

#ifndef TRJLIB_SAXSDEBYE_H_
#define TRJLIB_SAXSDEBYE_H_

#include "Saxs.h"
/** \brief Extends saxs to full Debey calculation through the ComputeSF method
 *
 */
class SaxsDebye: public Saxs {
public:
	SaxsDebye();
	SaxsDebye(double dq0, double qcut0): Saxs(dq0,qcut0){};
	SaxsDebye(int MyOrder,double dq0, double qcut0): Saxs(MyOrder,dq0,qcut0){};

	virtual void ComputeSAXS(RhoSaxs *,const MAtoms *);
	virtual ~SaxsDebye();
};

#endif /* TRJLIB_SAXSDEBEY_H_ */
