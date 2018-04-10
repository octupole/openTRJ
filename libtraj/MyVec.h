/*
 * MyVec.h
 *
 *  Created on: Jun 8, 2011
 *      Author: marchi
 */

#ifndef MYVEC_H_
#define MYVEC_H_
#include <iostream>
#include <Ftypedefs.h>

#define RAD2DEG		 (180.0/M_PI)

template <class T>
static inline T cos_angle(const T a[DIM],const T b[DIM])
{
  /* This version does not need the invsqrt lookup table */
  T   cosval;
  int    m;
  double aa,bb,ip,ipa,ipb; /* For accuracy these must be double! */

  ip=ipa=ipb=0.0;
  for(m=0; (m<DIM); m++) {		/* 18		*/
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cosval=ip/sqrt(ipa*ipb); 		/* 12		*/
					/* 30 TOTAL	*/
  if (cosval > 1.0)
    return  1.0;
  if (cosval <-1.0)
    return -1.0;

  return cosval;
}

template <class T> static inline T norm2(const T a[DIM])
{
  return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

template <class T> static inline T norm(const T a[DIM])
{
  return static_cast<T>(sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]));
}

template <class T>
struct myvector{
	T a[DIM];
	T & operator[](int i){
		try{
		if(i <DIM) return a[i];
		else throw " myvector boundaries reached ";
		}
		catch(const char * s){
			std::cout << s << std::endl;
			exit(1);
		}
	}
};

template <class T>
static inline myvector<T> Versor(const T a[DIM]){
	T mynorm=norm(a);
	myvector<T> temp;
	for(int i=0;i<DIM;i++) temp[i]=a[i]/mynorm;
	return temp;
}

#endif /* MYVEC_H_ */
