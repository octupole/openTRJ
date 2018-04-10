/*
 * DiffCoeffs.cpp
 *
 *  Created on: Jul 12, 2011
 *      Author: marchi
 */
#ifdef DIFFCOEFFS_H_

#include "DiffCoeffs.h"

template <unsigned int Order> int DiffCoeffs<Order>::nx=0;
template <unsigned int Order> int DiffCoeffs<Order>::ny=0;
template <unsigned int Order> int DiffCoeffs<Order>::nz=0;

template <unsigned int Order>
	DiffCoeffs<Order>::DiffCoeffs(int nx0, int ny0, int nz0, double dx, double dy, double dz)
		{nx=nx0;ny=ny0; nz=nz0;factx=float (nx0)/dx;facty=float (ny0)/dy;factz= float (nz0)/dz;}

template <unsigned int Order>
double DiffCoeffs<Order>::coeffs(int nn, int k){
	int ngrid[DIM]={nx,ny,nz};
	double grid[DIM]={factx,facty,factz};
	int n0=ngrid[nn];
	double fact;
	switch (Order){
	case 2:
		fact=coef[0]*sin(2.0*M_PI*(double) (k)/ (double) (n0));
		break;
	case 4:
		fact=coef[0]*sin(2.0*M_PI*(double) (k)/ (double) (n0))+coef[1]*sin(4.0*M_PI*(double) (k)/ (double) (n0));
		break;
	case 6:
		fact=coef[0]*sin(2.0*M_PI*(double) (k)/ (double) (n0))+coef[1]*sin(4.0*M_PI*(double) (k)/ (double) (n0))
			+coef[2]*sin(6.0*M_PI*(double) (k)/ (double) (n0));
		break;
	case 8:
		fact=coef[0]*sin(2.0*M_PI*(double) (k)/ (double) (n0))+coef[1]*sin(4.0*M_PI*(double) (k)/ (double) (n0))
			+coef[2]*sin(6.0*M_PI*(double) (k)/ (double) (n0))+coef[3]*sin(8.0*M_PI*(double) (k)/ (double) (n0));
		break;
	}
	return fact*grid[nn];
}



template <unsigned int Order> DiffCoeffs<Order>::DiffCoeffs() {

	// TODO Auto-generated constructor stub

}

template <unsigned int Order> DiffCoeffs<Order>::~DiffCoeffs() {
	// TODO Auto-generated destructor stub
}

#endif
