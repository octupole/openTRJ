/*
 * Saxs.h
 *
 *  Created on: Jun 23, 2015
 *      Author: marchi
 */

#ifndef SRC_SAXS_H_
#define SRC_SAXS_H_

#include "ScatteringFactors.h"
#include "ScatteringFactorsN.h"
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <functional>
#include <cmath>
#include "opgather.h"
#include "Atoms.h"

#include "histograms.hpp"
#include "integrate_function_adapt_simpson.h"
#include "RhoSaxs.h"
#include "RhoSaxsLI.h"
#include "specialfunctions.h"
#include "SaxsHistogram.h"
#include "SaxsHistogramSpline.h"
#include "myEnums.hpp"
#include "NewMPI.h"
#include "CenterMassBW3.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif
#include "stdafx.h"
#include "interpolation.h"
#include "specialfunctions.h"
#include "Spline1DInterpolant.h"



using namespace alglib;

using namespace Array;
using namespace fftwpp;
template <typename T>
inline char * as_byte(T & y){
	return reinterpret_cast<char *> (&y);
}
/** \brief Base class for the saxs calculation
 *         Methods for saxs computation and histogram build
 */
class Saxs {

protected:
	using AtomsD=Atoms<double>;
	using MAtoms=AtomsD;
	using MetricD=Metric<double>;
	using Dvect=DVECT::DDvect<double>;
	using Matrix=MATRIX::MMatrix<double>;
	using CenterMassD=CenterMass<double>;

	size_t align=sizeof(Complex);
	size_t count{0};
	double Vol{0.0};
	size_t Ntot{0},Ntot_in{0},Nsolute{0};
	double alpha{10.0};
	MetricD Mt;
	Matrix MCO{0.0};
	Matrix MOC{0.0};
	vector<std::function<double(double)> > MyFs;
	vector<size_t> MyNas;
	double dq=0.4,qcut=25.0;
	double dq_orig{0.0};
	double hx=0.05,mycut=4.0,mycut0=0.08;

	int order{-1};
	size_t nx{0},ny{0},nz{0},Nx{0},Ny{0},Nz{0};
	size_t nzp{0};
	SaxsHistogram Sgr, S, SContrast, Sdiff;
	map<const string,ScatteringFactors::opsfact> Sfacts;
	map<const string,vector<size_t> > iSfacts;
	int oneCount{1};
	vector<int>  Lst_i;
	vector<string> at1_i;
	void WriteIt(std::ostream &);
	array3<Complex> I_k;
	SaxsHistogram * qdfx{nullptr};
	SaxsHistogram * gofrx{nullptr};
	bool bAvg{false};
	bool first_time{true};
	bool noSplineOut{false};
	double SuperCell0{-1.0};
	double SuperCell{-1.0};
	double unitsR{10.0};
	double unitsQ{1.0/unitsR};
	double MassSolute{-1.0};
	bool bSans{false};
	void __copy(const Saxs &);
	virtual void Modulus(array3<Complex> &,array3<Complex> &);
	virtual SaxsHistogram __Qhistogram();
	virtual void __shift(AtomsD * y) {};
public:
	Saxs();
	Saxs(double dq0, double qcut0): dq{dq0}, dq_orig{dq0}, qcut{qcut0}{};
	Saxs(int MyOrder,double dq0, double qcut0): dq{dq0}, dq_orig{dq0},qcut{qcut0}, order{MyOrder}{};
	Saxs(const Saxs &);
	Saxs & operator-=(Saxs &);
	Saxs & operator=(const Saxs &);
	array3<Complex> getI_k(){return I_k;};
	void Clear();
	void SetupQdf();
	void SetMass(double a){MassSolute=a;};
	void SetSplineout(){noSplineOut=true;};
	bool bSplineout(){return noSplineOut;};
	void Averages();
	void Reduce(Parallel::NewMPI *);
	void Allocate(size_t,size_t,size_t);
	void CheckBin(double dq_min){if(dq_min>dq)dq=dq_min;}
	void setHist(double hhx, double cut){hx=hhx,mycut=cut;}
	void readContrast(std::ifstream &);
	size_t gNtot(){return Ntot;};
	size_t gNsolute(){return Nsolute;};
	Enums::Padding pickPadding(const vector<int> &) const;
	bool HaveSolute(vector<int> &) const ;
	bool HaveSuperCell(){return SuperCell0 > 0;}
	void bPrintIq(std::ofstream &);
	void bReadIq(std::ifstream &);
	void setSuperCell0(double x){SuperCell0=x;};
	void setSuperCell(double);
	void setbAvgTrue(){bAvg=true;}

	virtual void Setup(const vector<string> &,bool);
	virtual void Setup(const vector<int> &, const vector<string> &, bool);
	virtual void ComputeSq(RhoSaxs *,const MAtoms *);
	virtual void ComputeSAXS(RhoSaxs *,const MAtoms *);
	virtual void ComputeSANS(RhoSaxs *,const MAtoms *);
	void GofR();

	template <Compute myCompute>
	void Compute(RhoSaxs * ,MAtoms * );
	virtual ~Saxs();
	friend std::ostream & operator<<(std::ostream &, Saxs & );
};

#endif /* SRC_SAXS_H_ */

