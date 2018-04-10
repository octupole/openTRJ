/*
 * MyUtilClass.h
 *
 *  Created on: Jan 20, 2012
 *      Author: marchi
 */

#ifndef MYUTILCLASS_H_
#define MYUTILCLASS_H_

#include <iostream>
#include <iomanip>
#include <sstream>

#include <cmath>
#include <cstdlib>
#include <complex>

#include <Ftypedefs.h>

const double eps=1.0e-4;
const bool DEBUG=false;
using namespace Typedefs;

namespace MATRIX{
	template<class T>
	struct MMatrix;
}
namespace DVECT{
	template <class T>
	T sqr(T y){return y*y;};

	template <class T>
	struct DDvect{
		T x[DIM]{0,0,0};
		DDvect(){};
		DDvect(T x0, T y0, T z0){
			x[XX]=x0;
			x[YY]=y0;
			x[ZZ]=z0;
		}
		DDvect(const T y[DIM]){
			for(int o=0;o<DIM;o++) x[o]=y[o];
		}
		DDvect(T y){
			for(int o=0;o<DIM;o++) x[o]=y;
		}
		DDvect(const int y){
			for(int o=0;o<DIM;o++) x[o]=static_cast<T> (y);
		}
		T operator*(const DDvect & y){
			T temp=0.0;
			for(int o=0;o<DIM;o++) temp+=x[o]*y.x[o];
			return temp;
		};
		DDvect operator*(T y){
			DDvect temp;
			for(int o=0;o<DIM;o++){
				temp[o]=(*this)[o]*y;
			}
			return temp;
		};

		MATRIX::MMatrix<T> operator%(DDvect &);

		DDvect  operator^(const DDvect & y){
			DDvect temp;
			temp[XX]=x[YY]*y.x[ZZ]-x[ZZ]*y.x[YY];
			temp[YY]=x[ZZ]*y.x[XX]-x[XX]*y.x[ZZ];
			temp[ZZ]=x[XX]*y.x[YY]-x[YY]*y.x[XX];
			return temp;
		};
		DDvect  operator+(const DDvect & y){
			DDvect temp;
			for(int o=0;o<DIM;o++) temp[o]=x[o]+y.x[o];
			return temp;
		};
		DDvect  operator-(const DDvect & y){
			DDvect temp;
			for(int o=0;o<DIM;o++) temp[o]=x[o]-y.x[o];
			return temp;
		};
		DDvect  operator-(){
			DDvect temp;
			for(int o=0;o<DIM;o++) temp[o]=-x[o];
			return temp;
		};
		bool operator==(T y) const {
			bool is_y=true;
			for(int i=0;i<DIM;i++)
				if(x[i] != y) {is_y=false;break;}
			return is_y;
		}
		bool operator!=(T y) const {
			bool is_y=false;
			for(int i=0;i<DIM;i++)
				if(x[i] != y) {is_y=true;break;}
			return is_y;
		}
		DDvect & operator=(const DDvect & y){
			for(int o=0;o<DIM;o++) x[o]=y.x[o];
			return *this;
		};
		DDvect & operator=(T y){
			for(int o=0;o<DIM;o++) x[o]=y;
			return *this;
		};
		DDvect & operator=(const T y[DIM]){
			for(int o=0;o<DIM;o++) x[o]=y[o];
			return *this;
		};
		DDvect & operator()(const T y[DIM]){
			for(int o=0;o<DIM;o++) x[o]=y[o];
			return *this;
		};
		DDvect & operator()(const DDvect & y){
			*this=y;
			return *this;
		};
		DDvect & operator()(T x0, T y0, T z0){
			x[XX]=x0;
			x[YY]=y0;
			x[ZZ]=z0;
			return *this;
		};
		DDvect & operator+=(const DDvect & y){
			for(int o=0;o<DIM;o++) x[o]+=y.x[o];
			return *this;
		};
		DDvect & operator-=(const DDvect & y){
			for(int o=0;o<DIM;o++) x[o]-=y.x[o];
			return *this;
		};
		DDvect & operator/=(T y){
			for(int o=0;o<DIM;o++) x[o]/=y;
			return *this;
		}
		DDvect & operator*=(T y){
			for(int o=0;o<DIM;o++) x[o]*=y;
			return *this;
		}
		T Dist(DDvect & y){
			DDvect z=*this-y;
			return sqrt(z*z);
		}
		T Dist2(DDvect & y){
			DDvect z=*this-y;
			for(auto o=0;o<DIM;o++) z[o]=z[o]-rint(z[o]);
			return z*z;
		}
		void Displ(DDvect &, matrix &, matrix & );
		void Displ(DDvect &, MATRIX::MMatrix<T> & , MATRIX::MMatrix<T> & );
		DDvect Minus(DDvect &, MATRIX::MMatrix<T> & , MATRIX::MMatrix<T> & );
		DDvect Minus(DDvect &, matrix &, matrix & );
		void PBC(MATRIX::MMatrix<T> &, MATRIX::MMatrix<T> &);
		void PBC(matrix &, matrix &);
		T Dist(DDvect &, MATRIX::MMatrix<T> & , MATRIX::MMatrix<T> & );
		T Dist(DDvect &, matrix & , matrix & );

		DDvect & operator()(T y){
			for(int o=0;o<DIM;o++) x[o]=y;
			return *this;
		}
		DDvect operator/(T y){
			DDvect temp(*this);
			for(int o=0;o<DIM;o++) temp.x[o]/=y;
			return temp;
		}
		T & operator[](const int i){
			try {if(i>DIM-1) throw "Goes above dimensions set";}
			catch(const char * s){std::cout << s << std::endl;exit(1);}
			return x[i];
		};
		const T & operator[](const int i) const {
			try {if(i>DIM-1) throw "Goes above dimensions set";}
			catch(const char * s){std::cout << s << std::endl;exit(1);}
			return x[i];
		};
		void normalize(){
			T n=1.0/sqrt(x[XX]*x[XX]+x[YY]*x[YY]+x[ZZ]*x[ZZ]);
			*this*=n;
		}
		DDvect normal(){
			DDvect tmp=*this;
			T n=1.0/sqrt(x[XX]*x[XX]+x[YY]*x[YY]+x[ZZ]*x[ZZ]);
			tmp*=n;
			return tmp;
		}
		T Norm(){return sqrt(sqr(x[XX])+sqr(x[YY])+sqr(x[ZZ]));};
		T Norm2(){return sqr(x[XX])+sqr(x[YY])+sqr(x[ZZ]);};
		friend DDvect operator*(const T & y, const DDvect & z){
			DDvect temp;
			for(int o=0;o<DIM;o++)
					temp[o]=z[o]*y;
			return temp;
		};

		template<class T1>
		friend std::ostream & operator<<(std::ostream &, const DDvect<T1> &);

		template <class T1>
		friend std::stringstream & operator>>(std::stringstream & sin, DDvect<T1> &);

	};

}
namespace MATRIX{
	using namespace DVECT;

	template<class T>
	struct MMatrix{
		T m[DIM][DIM];
		MMatrix(){};
		MMatrix(T y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++) m[o][p]=y;
		}
		MMatrix(const int y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++) m[o][p]=static_cast<T> (y);
		}
		MMatrix(const T y[DIM][DIM]){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++) m[o][p]=y[o][p];
		};
		MMatrix(const DDvect<T> x, const DDvect<T> y, const DDvect<T> z){
			for(int o=0;o<DIM;o++) m[XX][o]=x[o];
			for(int o=0;o<DIM;o++) m[YY][o]=y[o];
			for(int o=0;o<DIM;o++) m[ZZ][o]=z[o];
		}
		T Trace() const {
			T tmp=0.0;
			for(int o=0;o<DIM;o++) tmp+=(*this)[o][o];
			return tmp;
		}
		MMatrix Transpose(){
			MMatrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=(*this)[p][o];
			return temp;
		}
		MMatrix Inversion(){
			MMatrix temp;
			T det;
			det=m[XX][XX]*(m[YY][YY]*m[ZZ][ZZ]-m[ZZ][YY]*m[YY][ZZ])-m[XX][YY]*(m[YY][XX]*m[ZZ][ZZ]-
					m[YY][ZZ]*m[ZZ][XX])+m[XX][ZZ]*(m[YY][XX]*m[ZZ][YY]-m[YY][YY]*m[ZZ][XX]);
			try{
				if(det == 0.0) throw std::string(" MMatrix is singular !!");
			}catch(const std::string & s){
				std::cout << s << std::endl;
				exit(-1);
			}
			temp[XX][XX]=(m[YY][YY]*m[ZZ][ZZ]-m[ZZ][YY]*m[YY][ZZ])/det;
			temp[YY][XX]=-(m[YY][XX]*m[ZZ][ZZ]-m[YY][ZZ]*m[ZZ][XX])/det;
			temp[ZZ][XX]=(m[YY][XX]*m[ZZ][YY]-m[ZZ][XX]*m[YY][YY])/det;
			temp[XX][YY]=-(m[XX][YY]*m[ZZ][ZZ]-m[XX][ZZ]*m[ZZ][YY])/det;
			temp[YY][YY]=(m[XX][XX]*m[ZZ][ZZ]-m[XX][ZZ]*m[ZZ][XX])/det;
			temp[ZZ][YY]=-(m[XX][XX]*m[ZZ][YY]-m[ZZ][XX]*m[XX][YY])/det;
			temp[XX][ZZ]=(m[XX][YY]*m[YY][ZZ]-m[XX][ZZ]*m[YY][YY])/det;
			temp[YY][ZZ]=-(m[XX][XX]*m[YY][ZZ]-m[YY][XX]*m[XX][ZZ])/det;
			temp[ZZ][ZZ]=(m[XX][XX]*m[YY][YY]-m[YY][XX]*m[XX][YY])/det;

			return temp;
		}
		const MMatrix & Unit(){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=(o==p)?1:0;
			return *this;
		};
		const T * operator[](const int i) const{
			try {if(i>DIM-1) throw "Goes above dimensions set";}
			catch(const char * s){std::cout << s << std::endl;exit(1);}
			return m[i];
		};
		T * operator[](const int i){
			try {if(i>DIM-1) throw "Goes above dimensions set";}
			catch(const char * s){std::cout << s << std::endl;exit(1);}
			return m[i];
		};
		MMatrix operator*(const MMatrix & y) const{
			MMatrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++){
					temp[o][p]=0.0;
					for(int r=0;r<DIM;r++)
						temp[o][p]+=(*this)[o][r]*y[r][p];
				}
			return temp;
		};
		DDvect<T> operator*(const DDvect<T> & y){
			DDvect<T> temp;
			for(int o=0;o<DIM;o++){
				temp[o]=0.0;
				for(int r=0;r<DIM;r++)
					temp[o]+=(*this)[o][r]*y[r];
			}
			return temp;
		};
		MMatrix RotTensor(const MMatrix & y){
			MMatrix & Rot=*this;
			MMatrix Rot_t=Rot.Transpose();
			MMatrix temp=y*Rot_t;
			MMatrix temp1=Rot*temp;
			return temp1;
		}
		const MMatrix operator+(const MMatrix & y){
			MMatrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=m[o][p]+y.m[o][p];
			return temp;
		};
		const MMatrix operator-(const MMatrix & y){
			MMatrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=m[o][p]-y.m[o][p];
			return temp;
		};
		const MMatrix operator-(){
			MMatrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=-m[o][p];
			return temp;
		};
		bool operator==(T y) const{
			bool is_y=true;
			for(int i=0;i<DIM;i++)
				for(int j=0;j<DIM;j++)
					if(m[i][j] != y) {is_y=false;break;}
			return is_y;
		}
		MMatrix & operator=(const MMatrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		MMatrix & operator=(MMatrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		MMatrix & operator=(T y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y;
			return *this;
		};
		MMatrix & operator=(const T  y[DIM][DIM]){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		MMatrix & operator()(const MMatrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		MMatrix & operator()(MMatrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]=y[o][p];
			return *this;
		};
		MMatrix & operator+=(const MMatrix & y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]+=y[o][p];
			return *this;
		};
		MMatrix operator/(T y){
			MMatrix temp(*this);
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp.m[o][p]/=y;
			return temp;
		}
		MMatrix & operator/=(T y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]/=y;
			return *this;
		};
		MMatrix & operator*=(T y){
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					m[o][p]*=y;
			return *this;
		};
		MMatrix operator*(T y){
			MMatrix temp;
			for(int o=0;o<DIM;o++)
				for(int p=0;p<DIM;p++)
					temp[o][p]=(*this)[o][p]*y;

			return temp;
		};

		template<class T1>
		friend MMatrix<T1> operator*(const T1 y, const MMatrix<T1> & z);

		template<class T1>
		friend void copyMatTomat(MMatrix<T1> & x,const MMatrix<T1> & y);

		template<class T1>
		friend std::ostream & operator<<(std::ostream &, const MMatrix<T1> &);
	};

}


#endif /* MYUTILCLASS_H_ */
