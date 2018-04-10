/*
 * Metric.h
 *
 *  Created on: May 25, 2011
 *      Author: marchi
 */

#ifndef METRIC_H_
#define METRIC_H_
#include <cmath>
#include <vector>

#include <Ftypedefs.h>
#include <MyUtilClass.h>
#include <MyVec.h>
using std::vector;

using namespace DVECT;
using namespace MATRIX;
template <typename T>
class Metric {
	using Matrix=MMatrix<T>;
	Matrix co, oc;
	const Matrix Idtity={{1,0,0},{0,1,0},{0,0,1}};
	void invertCO();

	T AA(){return norm(co[XX]);}
	T BB(){return norm(co[YY]);}
	T CC(){return norm(co[ZZ]);}
	T Alpha(){
		T alpha;
		if (norm2(co[YY])*norm2(co[ZZ])!=0)
			alpha = acos(cos_angle(co[YY],co[ZZ]));
		else
			alpha = 0.5*M_PI;
		return alpha;
	}
	T Beta(){
		T beta;
		if (norm2(co[XX])*norm2(co[ZZ])!=0)
			beta  = acos(cos_angle(co[XX],co[ZZ]));
		else
			beta  = 0.5*M_PI;
		return beta;
	}

	T Gamma(){
		T gamma;
		if (norm2(co[XX])*norm2(co[YY])!=0)
			gamma = acos(cos_angle(co[XX],co[YY]));
		else
			gamma = 0.5*M_PI;
		return gamma;
	}
public:
	Metric();
	Metric(const T,const T,const T,
			const T,const T,const T);
	Metric(const Matrix &);

	Metric(const Metric &);
	virtual ~Metric();
	Metric & operator()(const Matrix &);
	Metric & operator()(const Metric &);
	Metric & operator=(const Metric & );
	Metric & operator+=(const Metric &);
	Metric operator/(const T);
	Matrix getCO() {return co;};
	Matrix  getOC() {return oc;};

	T getVol(){ return (T) (co[XX][XX]*(co[YY][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[YY][ZZ])
		  -co[YY][XX]*(co[XX][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[XX][ZZ])
		  +co[ZZ][XX]*(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ]));};
	vector<T> getParas(){
		vector<T> d=vector<T>(DIM*2);
		d[0]=AA();
		d[1]=BB();
		d[2]=CC();
		d[3]=RAD2DEG*Alpha();
		d[4]=RAD2DEG*Beta();
		d[5]=RAD2DEG*Gamma();
		return d;}
};

#endif /* METRIC_H_ */
