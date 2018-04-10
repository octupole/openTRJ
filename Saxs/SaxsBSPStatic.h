/*
 * SaxsBSPStatic.h
 *
 *  Created on: May 10, 2016
 *      Author: marchi
 */

#ifndef SRC_SAXSBSPSTATIC_H_
#define SRC_SAXSBSPSTATIC_H_

#include "SaxsBSP.h"

class SaxsBSPStatic: public SaxsBSP {
protected:
	int time{0},nskip{1};
	vector<RhoSaxs> Rho_est;
public:
	SaxsBSPStatic();
	SaxsBSPStatic(const SaxsBSPStatic & y):SaxsBSP::SaxsBSP(y),time{y.time}{};
	SaxsBSPStatic(int mySkip,int MyOrder,double dq0, double qcut0): nskip{mySkip},SaxsBSP::SaxsBSP(MyOrder,dq0,qcut0){};

	virtual ~SaxsBSPStatic();
	void ComputeSAXS(RhoSaxs *,const MAtoms *);
	void ComputeSANS(RhoSaxs *,const MAtoms *);
};

#endif /* SRC_SAXSBSPSTATIC_H_ */
