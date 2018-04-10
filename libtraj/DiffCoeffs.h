/*
 * DiffCoeffs.h
 *
 *  Created on: Jul 12, 2011
 *      Author: marchi
 */

#ifndef DIFFCOEFFS_H_
#define DIFFCOEFFS_H_

#define DIM 3
#include <cmath>

template <unsigned int Order>
class DiffCoeffs {
	static int nx,ny,nz;
	double factx{0.0},facty{0.0},factz{0.0};
	static double coef[Order/2];
public:
	DiffCoeffs();
	DiffCoeffs(int nx0, int ny0, int nz0){nx=nx0;ny=ny0; nz=nz0;}
	DiffCoeffs(int nx0, int ny0, int nz0, double dx, double dy, double dz);
	void setD(double dx, double dy, double dz){factx=float (nx)/dx;facty=float (ny)/dy;factz= float (nz)/dz;};
	void setN(int nx0, int ny0, int nz0){nx=nx0;ny=ny0; nz=nz0;}
	double coeffs(int,int);
	virtual ~DiffCoeffs();
};

#include "DiffCoeffs.cpp"
#endif /* DIFFCOEFFS_H_ */
