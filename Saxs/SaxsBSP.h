/*
 * SaxsBSP.h
 *
 *  Created on: Jun 26, 2015
 *      Author: marchi
 */

#ifndef SRC_SAXSBSP_H_
#define SRC_SAXSBSP_H_

#include <functional>
#include "Saxs.h"
#include "BSpline.h"
#include "BSpmod.h"
#include "RhoSaxsBSP.h"

/** \brief Extends the saxs class for cardinal B-spline
 *
 */
class SaxsBSP: public Saxs {
protected:
	BSpmod * bsp_modx=nullptr;
	virtual void Modulus(array3<Complex> &,array3<Complex> &);
	virtual void __shift(AtomsD *);

public:
	SaxsBSP();
	SaxsBSP(int,double, double);

	virtual void Setup(const vector<string> &,bool);
	virtual void Setup(const vector<int> &, const vector<string> &,bool);

	virtual ~SaxsBSP();
	friend std::ofstream & operator<<(std::ofstream &, SaxsBSP & );
};

#endif /* SRC_SAXSBSP_H_ */
