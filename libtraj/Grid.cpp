/*
 * Grid.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: marchi
 */

#include "Grid.h"

template <unsigned int Ndim> const rvec Grid<Ndim>::cube[VERTEX]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};
template <unsigned int Ndim> int Grid<Ndim>::SetNo=0;
static int Fcounter=0;

template <class T>
T sqr(T a){return a*a;}

template <unsigned int Ndim>
Grid<Ndim>::Grid()
{
		try{
			if(nnx+nny+nnz==0) throw "\nWarning: Metrics initialization is postponed\n";
		}
		catch(const char * s){
			std::cout << s << std::endl;
		}

}
template <unsigned int Ndim>
Grid<Ndim>::Grid(size_t nx,size_t ny, size_t nz){
	Allocate(nx,ny,nz);
}
template <unsigned int Ndim>
Grid<Ndim>::Grid(Grid<Ndim> & y) {
	nnx=y.nnx;nny=y.nny;nnz=y.nnz;
	Allocate();
	array4<double>::operator=(y);
}
template <unsigned int Ndim>
Grid<Ndim>::Grid(const Grid<Ndim> & y){
	nnx=y.nnx;nny=y.nny;nnz=y.nnz;
	Allocate();
	array4<double>::operator=(y);
}
template <unsigned int Ndim>
Grid<Ndim>::Grid(const double & x){
	this->Allocate();
	*this=x;
}

template <unsigned int Ndim>
inline Grid<Ndim> & Grid<Ndim>::operator=(const double y){
	double * pt=&(*this)[0][0][0][0];
	for(int i=0;i<static_cast<int>(Ndim*nnx*nny*nnz);i++) pt[i]=y;
	return *this;
}

template <unsigned int Ndim>
inline Grid<Ndim> & Grid<Ndim>::operator=(const Grid<Ndim> & y){
	nnx=y.nnx;nny=y.nny;nnz=y.nnz;
	Allocate();
	array4<double>::operator=(y);
	return *this;
}

template <unsigned int Ndim>
inline Grid<Ndim> & Grid<Ndim>::operator+=(const Grid<Ndim> & y){
	double * pt=&(*this)[0][0][0][0];
	double * py=&y[0][0][0][0];
	for(int i=0;i<static_cast<int>(Ndim*nnx*nny*nnz);i++) pt[i]+=py[i];
	return *this;
}

template <unsigned int Ndim>
inline Grid<Ndim> Grid<Ndim>::operator-(Grid<Ndim> & a){
	Grid<Ndim> temp=*this;
	temp-=a;
	return temp;
}

template <unsigned int Ndim>
inline Grid<Ndim> Grid<Ndim>::operator+(Grid<Ndim> & a){
	Grid<Ndim> temp=*this;
	temp+=a;
	return temp;
}

template <unsigned int Ndim>
inline Grid<Ndim> Grid<Ndim>::operator*(const double & a){
	Grid<Ndim> temp=*this;
	temp*=a;
	return temp;
}

template <unsigned int Ndim>
inline Grid<Ndim> Grid<Ndim>::operator/(const double & a){
	Grid<Ndim> temp=*this;
	temp/=a;
	return temp;
}
template <unsigned int Ndim>
inline Grid<Ndim> Grid<Ndim>::operator/(Grid<Ndim> & a){
	Grid<Ndim> temp=*this;
	for(unsigned int l=0;l<Ndim;l++)
		for(unsigned int i=0;i<nnx;i++)
			for(unsigned int j=0;j<nny;j++)
				for(unsigned int k=0;k<nnz;k++){
					if(a[l][i][j][k] == 0) temp[l][i][j][k]=0.0;
					else temp[l][i][j][k]/=a[l][i][j][k];
				}
	return temp;
}
template <unsigned int Ndim>
inline void Grid<Ndim>::setMetric(const Matrix & co_in){
	matrix co;
	for(int n=0;n<DIM;n++)
		for(int m=0;m<DIM;m++)
			co[n][m]=co_in[n][m];
	setMetric(co);
}
template <unsigned int Ndim>
inline void Grid<Ndim>::setMetric(const matrix & co_in)
{
  /* Save some time by assuming lower right part is zero */
  for(int i=0;i<DIM;i++)
	  for(int j=0;j<DIM;j++) co[i][j]=co_in[i][j];

  float tmp=1.0/(co[XX][XX]*co[YY][YY]*co[ZZ][ZZ]);
  oc[XX][XX]=co[YY][YY]*co[ZZ][ZZ]*tmp;
  oc[YY][XX]=0;
  oc[ZZ][XX]=0;
  oc[XX][YY]=-co[XX][YY]*co[ZZ][ZZ]*tmp;
  oc[YY][YY]=co[XX][XX]*co[ZZ][ZZ]*tmp;
  oc[ZZ][YY]=0;
  oc[XX][ZZ]=(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ])*tmp;
  oc[YY][ZZ]=-co[YY][ZZ]*co[XX][XX]*tmp;
  oc[ZZ][ZZ]=co[XX][XX]*co[YY][YY]*tmp;
  Volume=(co[XX][XX]*(co[YY][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[YY][ZZ])
		  -co[YY][XX]*(co[XX][YY]*co[ZZ][ZZ]-co[ZZ][YY]*co[XX][ZZ])
		  +co[ZZ][XX]*(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ]));
}


template <unsigned int Ndim>
inline double Grid<Ndim>::Alpha(){
	double alpha;
	if (norm2(co[YY])*norm2(co[ZZ])!=0)
		alpha = acos(cos_angle(co[YY],co[ZZ]));
	else
		alpha = 0.5*M_PI;

	return alpha;
}

template <unsigned int Ndim>
inline double Grid<Ndim>::Beta(){
	double beta;
	if (norm2(co[XX])*norm2(co[ZZ])!=0)
		beta  = acos(cos_angle(co[XX],co[ZZ]));
	else
		beta  = 0.5*M_PI;
	return beta;
}

template <unsigned int Ndim>
inline double Grid<Ndim>::Gamma(){
	double gamma;
	if (norm2(co[XX])*norm2(co[YY])!=0)
		gamma = acos(cos_angle(co[XX],co[YY]));
	else
		gamma = 0.5*M_PI;
	return gamma;
}

template <unsigned int Ndim>
Grid<Ndim>::~Grid() {
	// TODO Auto-generated destructor stub
}
template <unsigned int Ndim>
Grid<Ndim>::Grid(array4<Complex> & x){
	Allocate();
	for(unsigned int l=0;l<Ndim;l++)
		for(unsigned int i=0;i<nnx;i++)
			for(unsigned int j=0;j<nny;j++)
				for(unsigned int k=0;k<nnz;k++){
					(*this)[l][i][j][k]=x[l][i][j][k].real();
				}
}
template <unsigned int Ndim>
Grid<Ndim> & Grid<Ndim>::operator=(const array4<Complex> & a){
	try{
		if(this->Size() != a.Size()) throw "Grid and array3 arrays do not have the same size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(unsigned l=0;l<Ndim; l++){
		for(unsigned i=0;i<nnx;i++)
			for(unsigned j=0;j<nny;j++)
				for(unsigned k=0;k<nnz;k++)
					(*this)[l][i][j][k]=a[l][i][j][k].real();
				}
	return *this;

}
template <unsigned int Ndim>
Grid<Ndim>::Grid(array4<double> & x){
	Allocate();
	for(unsigned int l=0;l<Ndim;l++)
		for(unsigned int i=0;i<nnx;i++)
			for(unsigned int j=0;j<nny;j++)
				for(unsigned int k=0;k<nnz;k++){
					(*this)[l][i][j][k]=x[l][i][j][k];
				}
}
template <unsigned int Ndim>
inline Grid<Ndim> & Grid<Ndim>::operator=(const array4<double> & a){
	try{
		if(this->Size() != a.Size()) throw "Grid and array3 arrays do not have the same size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(unsigned l=0;l<Ndim; l++){
		for(unsigned i=0;i<nnx;i++)
			for(unsigned j=0;j<nny;j++)
				for(unsigned k=0;k<nnz;k++)
					(*this)[l][i][j][k]=a[l][i][j][k];
				}
	return *this;
}

template <unsigned int Ndim>
inline void Grid<Ndim>::TriLinear(Real gx, Real gy, Real gz, Real q[Ndim]){
	int n;
	for(n=0;n<VERTEX;n++){
		Real dd;
		dd=TriDiff(gx,cube[n][XX])*TriDiff(gy,cube[n][YY])*TriDiff(gz,cube[n][ZZ]);
		for(unsigned int o=0;o<Ndim;o++)
			qq[n][o]=dd*q[o];
	}
}

template <unsigned int Ndim>
inline Grid<Ndim>  operator-(const Grid<Ndim> & y){
	Grid<Ndim> temp;
	for(unsigned l=0;l<Ndim; l++){
		for(unsigned i=0;i<y.nnx;i++)
			for(unsigned j=0;j<y.nny;j++)
				for(unsigned k=0;k<y.nnz;k++)
					temp[l][i][j][k]=-y[l][i][j][k];
				}
	return temp;
}
template <class T, unsigned int Ndim>
inline Grid<Ndim> operator*(const T & a,const Grid<Ndim> & y){
	Grid<Ndim> temp;
	for(unsigned l=0;l<Ndim; l++){
		for(unsigned i=0;i<y.nnx;i++)
			for(unsigned j=0;j<y.nny;j++)
				for(unsigned k=0;k<y.nnz;k++)
					temp[l][i][j][k]=static_cast<double> (a)*y[l][i][j][k];
				}
	return temp;
}


template <>
inline Grid<1> & Grid<1>::operator=(const array3<double> & a){
	try{
		if(this->Size() != a.Size()) throw "Grid and array3 arrays do not have the same size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(unsigned int i=0;i<nnx;i++)
		(*this)[0][i]=a[i];
	return *this;
}
template <>
inline Grid<1>::Grid(array3<Complex> & x){
	Allocate();
	for(unsigned int i=0;i<nnx;i++)
		for(unsigned int j=0;j<nny;j++)
			for(unsigned int k=0;k<nnz;k++)
				(*this)[0][i][j][k]=x[i][j][k].real();
}

template <>
inline Grid<1>::Grid(array3<double> & x){
	Allocate();
	(*this)[0]=x;
}


template <>
inline Grid<1> & Grid<1>::operator=(const array3<Complex> & a){
	try{
		if(this->Size() != a.Size()) throw "Grid and array3 arrays do not have the same size ";
	}

	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(unsigned i=0;i<nnx;i++)
		for(unsigned j=0;j<nny;j++)
			for(unsigned k=0;k<nnz;k++){
				(*this)[0][i][j][k]=a[i][j][k].real();
			}
	return *this;
}

template class Grid<1>;
template class Grid<DIM>;
