/*
 * RhoSaxs.h
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#ifndef SRC_RHOSAXS_H_
#define SRC_RHOSAXS_H_

#include <map>
#include "BSpline.h"
#include "LagrangeInterpolation.h"
#include "VecRotate.h"
#include "myEnums.hpp"
#include "Grid.h"
#include "Atoms.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "Finalize.h"

using std::map;
const int NMult=1;

using namespace Enums;
/** \brief Derived class of Grid<1>
 *         A tridimensional grid storing the atomic density.
 */
class RhoSaxs: public Grid<1> {
protected:
	using Matrix=MMatrix<double>;
	using Dvect=DDvect<double>;
	using MetricD=Metric<double>;
	using AtomsD=Atoms<double>;
	using MAtoms=AtomsD;
	const double SHELL{0.07};
	static bool firstTime;
	map<string,double> MyPadding;
	Padding exPadding{Enums::avgDensity};
	string Type{" "};
	void zeroPadding(RhoSaxs *,int,int,int);
	void avgPadding(RhoSaxs *,int,int,int);
	void myPadding(RhoSaxs *,int,int,int);
	void PBCPadding(RhoSaxs *,int,int,int);
	virtual void __Density(const int x,const AtomsD * y, vector<size_t> & z, string w
			, vector<double> & wei);
public:
	RhoSaxs(){};
	RhoSaxs(const RhoSaxs & y): Grid<1>::Grid<1>(y){}
	RhoSaxs(size_t nx,size_t ny,size_t nz):Grid<1>::Grid<1>(nx,ny,nz) {};
	void MakeAvg();
	virtual RhoSaxs & operator=(const double y){
		this->Grid<1>::operator =(y);
		return *this;
	};
	virtual RhoSaxs & operator=(const RhoSaxs & y){
		this->Grid<1>::operator =(y);
		return *this;
	};
	void setPadding(const map<string,double> &);
	void selectPadding(Padding );

	/** \brief Function template doing padding for the supercell
	 *
	 * @param y Pointer to an instantiated RhoSaxs grid
	 * @param Nx Grid dimension on x
	 * @param Ny Grid dimension on y
	 * @param Nz Grid dimension on z
	 */
	template <Padding myPadd>
	void doPadding(RhoSaxs * y,int Nx, int Ny, int Nz){
		switch(myPadd){
		case zero:
			zeroPadding(y,Nx,Ny,Nz);
			break;
		case avgDensity:
			avgPadding(y,Nx,Ny,Nz);
			break;
		case myDensity:
			myPadding(y,Nx,Ny,Nz);
			break;
		case Periodic:
			PBCPadding(y,Nx,Ny,Nz);
			break;
		default:
			break;
		}
	}
//	void Extend(RhoSaxs *,int,Padding);
	/**
	 *
	 * @param x The order of the interpolant. Unused here.
	 * @param y Pointer to an instantiation of Atoms or its derived classes
	 * @param z Indirect addressing to specific particles meshed by the class
	 * @param w Unused
	 */
	void Density(const int x,const AtomsD * y, vector<size_t> & z, string w);
	void Density(const int x,const AtomsD * y, vector<size_t> & z, string w, vector<double> & wei);
	map<string,double> & getPadding(){return MyPadding;}
//	Enums::Padding WhichPadding(){if(MyPadding.empty()) return avgDensity; return myDensity;}
	Enums::Padding WhichPadding(){return exPadding;};
	virtual ~RhoSaxs();
};

#endif /* SRC_RHOSAXS_H_ */
