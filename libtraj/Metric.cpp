/*
 * Metric.cpp
 *
 *  Created on: May 25, 2011
 *      Author: marchi
 */

#include "../libtraj/Metric.h"
template <typename T>
Metric<T>::Metric(){
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]=Idtity[i][j];
			oc[i][j]=Idtity[i][j];
		}

}
template <typename T>
Metric<T>::Metric(const T a,const T b,const T c,
		const T alfa,const T beta,const T gamma){
	T degrad=M_PI/180.0;

	T ax=a;
	T alf=cos(degrad*alfa);
	T bet=cos(degrad*beta);
	T qt=sin(degrad*gamma);
	T gam=cos(degrad*gamma);
	T bx=b*gam;
	T by=b*qt;
	T cx=c*bet;
	T cy=c*(alf-bet*gam)/qt;
	T cz=sqrt(c*c-cx*cx-cy*cy);
	co[YY][XX]=0.0;
	co[ZZ][XX]=0.0;
	co[ZZ][YY]=0.0;
	co[XX][XX]=ax;
	co[XX][YY]=bx;
	co[XX][ZZ]=cx;
	co[YY][YY]=by;
	co[YY][ZZ]=cy;
	co[ZZ][ZZ]=cz;
}

template <typename T>
Metric<T>::Metric(const Matrix & co_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) co[i][j]=co_in[i][j];
	invertCO();
}
template <typename T>
Metric<T>::Metric(const Metric & Mt_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]=Mt_in.co[i][j];
			oc[i][j]=Mt_in.oc[i][j];
		}

}
template <typename T>
Metric<T> Metric<T>::operator/(const T a){
	Metric tmp=*this;
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			tmp.co[i][j]=tmp.co[i][j]/(a);
			tmp.oc[i][j]=tmp.oc[i][j]/(a);
		}
	return tmp;
}
template <typename T>
Metric<T> & Metric<T>::operator()(const Metric & Mt_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]=Mt_in.co[i][j];
			oc[i][j]=Mt_in.oc[i][j];
		}
	return *this;
}
template <typename T>
Metric<T> & Metric<T>::operator=(const Metric & Mt_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]=Mt_in.co[i][j];
			oc[i][j]=Mt_in.oc[i][j];
		}
	return *this;
}
template <typename T>
Metric<T> & Metric<T>::operator+=(const Metric & Mt_in) {
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) {
			co[i][j]+=Mt_in.co[i][j];
			oc[i][j]+=Mt_in.oc[i][j];
		}
	return *this;
}
template <typename T>
Metric<T> & Metric<T>::operator()(const Matrix & co_in){
	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++) co[i][j]=co_in[i][j];
	invertCO();
	return *this;
}
template <typename T>
Metric<T>::~Metric() {

}

template <typename T>
void Metric<T>::invertCO()
{
  /* Save some time by assuming lower right part is zero */
  real tmp=1.0/(co[XX][XX]*co[YY][YY]*co[ZZ][ZZ]);
  oc[XX][XX]=co[YY][YY]*co[ZZ][ZZ]*tmp;
  oc[YY][XX]=0;
  oc[ZZ][XX]=0;
  oc[XX][YY]=-co[XX][YY]*co[ZZ][ZZ]*tmp;
  oc[YY][YY]=co[XX][XX]*co[ZZ][ZZ]*tmp;
  oc[ZZ][YY]=0;
  oc[XX][ZZ]=(co[XX][YY]*co[YY][ZZ]-co[YY][YY]*co[XX][ZZ])*tmp;
  oc[YY][ZZ]=-co[YY][ZZ]*co[XX][XX]*tmp;
  oc[ZZ][ZZ]=co[XX][XX]*co[YY][YY]*tmp;
}
template class Metric<float>;
template class Metric<double>;
