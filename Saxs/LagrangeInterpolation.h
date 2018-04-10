/*
 * LagrangeInterpolation.h
 *
 *  Created on: Jul 7, 2015
 *      Author: marchi
 */

#ifndef SRC_LAGRANGEINTERPOLATION_H_
#define SRC_LAGRANGEINTERPOLATION_H_
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <functional>
#include <sstream>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::map;
using std::function;
using std::stringstream;



struct opLagrange{

};
/** \brief Lagrange interpolation class.
 * 		   Return the interpolant as a map of lambda functions
 */
class LagrangeInterpolation {
	static map<int,map<int,function<double(double)> > > Interpolant;
	map<int,function<double(double)> > MyPoly;
	size_t order=2;
	const size_t MaxOrder{6};
	vector<double> x;
public:
	LagrangeInterpolation();
	LagrangeInterpolation(int);
	const map<int,function<double(double)> > & Poly() const {return MyPoly;}
	function<double(double)> & operator[](int n){return MyPoly[n];}
	virtual ~LagrangeInterpolation();
};

#endif /* SRC_LAGRANGEINTERPOLATION_H_ */
