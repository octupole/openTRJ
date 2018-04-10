/*
 * BSpmod.h
 *
 *  Created on: Jun 26, 2015
 *      Author: marchi
 */

#ifndef SRC_BSPMOD_H_
#define SRC_BSPMOD_H_
#include "Ftypedefs.h"
#include <vector>
#include <string>
#include "BSpline.h"
#include <iostream>
#include <cmath>


using std::vector;
using std::string;
using std::cout;
using std::endl;
const auto KCUT=50;

struct bspmoduli{
	vector<double> x;
	vector<double> y;
	vector<double> z;
	bspmoduli & operator()(size_t nx,size_t ny,size_t nz){
		x=vector<double>(nx);
		y=vector<double>(ny);
		z=vector<double>(nz);
		return *this;
	}
};

class BSpmod {
	bspmoduli BSp;
	static size_t nx,ny,nz,order;
	static size_t ndim;
	void load_moduli();
	vector<double> DFTmod(const vector<double> &, size_t);
	void Gamma(vector<double> &);
	void Inverse();
public:
	BSpmod(size_t nnx,size_t nny,size_t nnz){
		nx=nnx;ny=nny;nz=nnz;
		auto Max3=[](size_t a,size_t b, size_t c)->size_t{size_t max=a<b?b:a;return max=max<c?c:max;};
		ndim=Max3(nx,ny,nz);
		BSp(nx,ny,nz);
		order=BSpline::Order();
		try{
			if(order == 0) throw string("Cannot run with a spline order equal to zero.");
		} catch(const string & s){
			cout << s <<endl;
			exit(0);
		}
		load_moduli();
	}
	const double & ModuliX(size_t o) const {return BSp.x[o];};
	const double & ModuliY(size_t o) const {return BSp.y[o];};
	const double & ModuliZ(size_t o) const {return BSp.z[o];};
	virtual ~BSpmod();
};

#endif /* SRC_BSPMOD_H_ */
