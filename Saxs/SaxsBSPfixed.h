/*
 * SaxsBSPfixed.h
 *
 *  Created on: Feb 8, 2018
 *      Author: marchi
 */

#ifndef SRC_SAXSBSPFIXED_H_
#define SRC_SAXSBSPFIXED_H_

#include "SaxsBSP.h"

class SaxsBSPfixed: public SaxsBSP {
	using SaxsBSP::SaxsBSP;
	int time{0},nskip{1};
	vector<RhoSaxs> Rho_est;
	RhoSaxs * Rho_fake;
public:
	virtual void ComputeSAXS(RhoSaxs * ,const MAtoms *);
	virtual ~SaxsBSPfixed();
};

#endif /* SRC_SAXSBSPFIXED_H_ */
