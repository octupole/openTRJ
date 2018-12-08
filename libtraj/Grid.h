/*
 * Grid.h
 *
 *  Created on: Jun 6, 2011
 *      Author: marchi
 */

#ifndef GRID_H_
#define GRID_H_
//#include "vec.h"
#include "Array.h"
#include <cstdlib>
#include <cmath>
#include "MyVec.h"
#include <vector>
#include <string>
#include <sstream>
#include "MyUtilClass.h"
#include <complex>
using Complex=std::complex<double>;


typedef float Real;
using namespace Array;

using std::vector;
using std::string;
using std::cout;
using std::endl;
using Typedefs::matrix;
using Typedefs::rvec;

#define VERTEX 8

template <unsigned int Ndim>
class Grid: public array4<double> {
protected:
	using Matrix=MATRIX::MMatrix<double>;
	using Dvect=DVECT::DDvect<double>;

	std::string name;
	matrix co{{0,0,0},{0,0,0},{0,0,0}},oc{{0,0,0},{0,0,0},{0,0,0}};
	double Volume{0.0};
	size_t nnx{0},nny{0},nnz{0};
	Real qq[VERTEX][Ndim];
	static const rvec cube[VERTEX];
	void setMetric(const matrix &);
	void setMetric(const Matrix &);
	double Alpha();
	double Beta();
	double Gamma();
	Real TriDiff(Real g, Real cubi){return((g-cubi>=0.0)? cubi-g+1 : g-cubi+1);};
	void TriLinear(Real, Real, Real, Real [Ndim]);
public:
	static int SetNo;
	Grid();
	Grid(Grid &);
	Grid(const Grid &);
	Grid(array3<Complex> &);
	Grid(array4<Complex> &);
	Grid(array4<double> &);
	Grid(array3<double> &);
	Grid(const double &);
	Grid(size_t,size_t,size_t);
	Grid & operator=(const double);
	Grid & operator=(const Grid &);
	Grid & operator=(const array3<double> &);
	Grid & operator=(const array3<Complex> &);
	Grid & operator=(const array4<Complex> &);
	Grid & operator=(const array4<double> &);
	Grid & operator+=(const Grid &);
	Grid operator*(const double &);
	Grid operator/(const double &);
	Grid operator/(Grid &);
	Grid operator-(Grid &);
	Grid operator+(Grid &);
	void setNN(unsigned int nx0, unsigned int ny0,unsigned int nz0){nnx=nx0;nny=ny0;nnz=nz0;};
	void set(const matrix & CO,unsigned int nx0=0,unsigned int ny0=0,unsigned int nz0=0){
		if(nx0+ny0+nz0) setNN(nx0,ny0,nz0);setMetric(CO);
	};

	void Allocate(){
		try{
			if(nnx==0 && nny==0&&nnz==0) throw "Cannot allocate if a dimension is zero!!";
		}
		catch(const char * s){
			std::cout << s << std::endl;
			exit(1);
		}
		if(!this->allocated) array4<double>::Allocate(Ndim,nnx,nny,nnz);
		else {
			this->Deallocate();
			array4<double>::Allocate(Ndim,nnx,nny,nnz);
		}
		*this=0.0;
	};
	void Allocate(size_t nx,size_t ny, size_t nz){
		nnx=nx;nny=ny;nnz=nz;
		if(!this->allocated) array4<double>::Allocate(Ndim,nnx,nny,nnz);
		else {
			this->Deallocate();
			array4<double>::Allocate(Ndim,nnx,nny,nnz);
		}
		*this=0.0;
	};
	double getVol(){return Volume;};
	matrix & getCO(){return co;};
	matrix & getOC(){return oc;};
	double getDV(){return Volume/static_cast<double>(nnx*nny*nnz);};
	unsigned int getnnx(){return nnx;};
	unsigned int getnny(){return nny;};
	unsigned int getnnz(){return nnz;};
	unsigned int getn0(){return Ndim;};
	unsigned int getnx(){return ny;};
	unsigned int getny(){return nz;};
	unsigned int getnz(){return nw;};
	void Xplor(FILE *,double=1.0);
	vector<double> Rdf(FILE * fp,const double [DIM], std::string, const double & = 2.0, const double & = 0.02);
	void RdfK(FILE * fp, std::string, const double & = 200.0, const double & = 0.02);

	void setname(std::string & nme){name=nme;};
	template <unsigned int Ndima>
	friend Grid<Ndima>  operator-(const Grid<Ndima> &);
	template <class T,unsigned int Ndima>
	friend Grid<Ndima>  operator*(const T &,const Grid<Ndima> &);
	virtual ~Grid();
};

#endif /* GRID_H_ */
