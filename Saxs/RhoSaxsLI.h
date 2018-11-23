/*
 * RhoSaxsLI.h
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#ifndef SRC_RHOSAXSLI_H_
#define SRC_RHOSAXSLI_H_

#include "RhoSaxs.h"

/** \brief Extends RhoSaxs class to do Lagrangian interpolation
 *
 */
class RhoSaxsLI: public RhoSaxs {
	virtual void __Density(const int x,const AtomsD * y, vector<size_t> & z, string w
			, vector<double> & wei);
public:
	RhoSaxsLI();
	virtual ~RhoSaxsLI();
};

#endif /* SRC_RHOSAXSLI_H_ */
