/*
 * RhoSaxsBSP.h
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#ifndef SRC_RHOSAXSBSP_H_
#define SRC_RHOSAXSBSP_H_

#include "RhoSaxs.h"
/** \brief Extends RhoSaxs class to do cardinal B-spline
 *
 */
class RhoSaxsBSP: public RhoSaxs {
	void __Density(const int x,const AtomsD * y, vector<size_t> & z, string w
			, vector<double> & wei);
public:
	RhoSaxsBSP(){};
	RhoSaxsBSP(const RhoSaxsBSP & y): RhoSaxs::RhoSaxs(y){};
	virtual ~RhoSaxsBSP();
};

#endif /* SRC_RHOSAXSBSP_H_ */
