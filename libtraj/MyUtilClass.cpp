/*
 * MyUtilClass.cpp
 *
 *  Created on: Jan 20, 2012
 *      Author: marchi
 */


#include "../libtraj/MyUtilClass.h"
namespace DVECT{

using namespace MATRIX;

	template<class T>
	MMatrix<T> DDvect<T>::operator %(DDvect<T> & y){
		MMatrix<T> temp;
		for(int o=0;o<DIM;o++)
			for(int p=0;p<DIM;p++)
				temp[o][p]=this->x[o]*y[p];
		return temp;
	};
	template<class T>
	std::ostream & operator<<(std::ostream & out, const DDvect<T> & z){
		out << std::setw(11) << std::setprecision(4) << std::fixed << z[XX] << " ";
		out << std::setw(11) << std::setprecision(4) << std::fixed << z[YY] << " ";
		out << std::setw(11) << std::setprecision(4) << std::fixed << z[ZZ] << std::endl;;

		return out;
	}
	template<class T>
	std::stringstream & operator>>(std::stringstream & sin, DDvect<T> & z){
		sin>>z[XX]>>z[YY]>>z[ZZ];
		return sin;
	}

	template<class T>
	MMatrix<T> operator%(const DDvect<T> & x,const DDvect<T> & y){
		MMatrix<T> temp;
		for(int o=0;o<DIM;o++)
			for(int p=0;p<DIM;p++)
				temp[o][p]=x[o]*y[p];
		return temp;
	};
	template<class T>
	T DDvect<T>::Dist(DDvect<T> & y, MMatrix<T> & co, MMatrix<T> & oc){
		DDvect<T> xc=Minus(y,co,oc);
		return sqrt(xc*xc);
	}

	template<class T>
	DDvect<T> DDvect<T>::Minus(DDvect<T> & y, MMatrix<T> & co, MMatrix<T> & oc){
		DDvect xa=oc*x-oc*y;
		DDvect xb;
		xb[XX]=xa[XX]-rint(xa[XX]);
		xb[YY]=xa[YY]-rint(xa[YY]);
		xb[ZZ]=xa[ZZ]-rint(xa[ZZ]);
		return co*xb;
	}
	template<class T>
	void DDvect<T>::Displ(DDvect<T> & y,MMatrix<T> & co, MMatrix<T> & oc){
		DDvect<T> xa=oc*x-oc*y;
		DDvect<T> xb=oc*x;
		xb[XX]-=rint(xa[XX]);
		xb[YY]-=rint(xa[YY]);
		xb[ZZ]-=rint(xa[ZZ]);
		*this=co*xb;
	}
	template<class T>
	void DDvect<T>::PBC(MMatrix<T> & co, MMatrix<T> & oc){
		DDvect<T> xa=oc*(*this);
		DDvect<T> xb;
		xb[XX]=xa[XX]-rint(xa[XX]);
		xb[YY]=xa[YY]-rint(xa[YY]);
		xb[ZZ]=xa[ZZ]-rint(xa[ZZ]);
		*this=co*xb;
	}
	template struct DDvect<float>;
	template struct DDvect<double>;
	template std::ostream & operator<<(std::ostream &, const DDvect<double> &);
	template std::ostream & operator<<(std::ostream &, const DDvect<float> &);
	template std::stringstream & operator>>(std::stringstream & , DDvect<float> &);
	template std::stringstream & operator>>(std::stringstream & , DDvect<double> &);

}
namespace MATRIX{
template<class T>
void copyMatTomat(MMatrix<T> & x,const MMatrix<T> & y){
	for(auto o=0;o<DIM;o++)
		for(auto p=0;p<DIM;p++)
			x[o][p]=y[o][p];
}


template<class T>
MMatrix<T> operator*(const T y, const MMatrix<T> & z){
	MMatrix<T> temp;
	for(int o=0;o<DIM;o++)
		for(int p=0;p<DIM;p++)
			temp[o][p]=z[o][p]*y;
	return temp;
};
template<class T>
	std::ostream & operator<<(std::ostream & out, const MMatrix<T> & z){
		out << "      MMatrix        " << std::endl;
		out << std::setw(11) << std::setprecision(4) << std::fixed << z[XX][XX] << " " << z[XX][YY] << " " << z[XX][ZZ] << std::endl;
		out << std::setw(11) << std::setprecision(4) << std::fixed<< z[YY][XX] << " " << z[YY][YY] << " " << z[YY][ZZ] << std::endl;
		out << std::setw(11) << std::setprecision(4) << std::fixed<< z[ZZ][XX] << " " << z[ZZ][YY] << " " << z[ZZ][ZZ] << std::endl;
		return out;
	}

template struct MMatrix<float>;
template struct MMatrix<double>;
template MMatrix<float> operator*(const float, const MMatrix<float> &);
template MMatrix<double> operator*(const double, const MMatrix<double> &);
template std::ostream & operator<<(std::ostream &, const MMatrix<float> &);
template std::ostream & operator<<(std::ostream &, const MMatrix<double> &);
template void copyMatTomat(MMatrix<float> & x,const MMatrix<float> & y);
template void copyMatTomat(MMatrix<double> & x,const MMatrix<double> & y);

}

