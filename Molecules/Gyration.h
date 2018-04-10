/*
 * Gyration.h
 *
 *  Created on: Aug 11, 2013
 *      Author: marchi
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include "json.hpp"
#include <functional>


#include "MyUtilClass.h"

#ifndef GYRATION_H_
#define GYRATION_H_
using namespace std;
using namespace DVECT;
using namespace MATRIX;
using nlohmann::json;
template <typename T>
class Gyration{
protected:
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	static json myJson;
	std::function<T(T)> mySqrt=[](T x){return x>0?sqrt(x):-1;};
	map<string,size_t> str_hash;
	double Radg{0};
	Dvect I,G,axis;
	Dvect * xcmx{nullptr};
	static double time_c;
	virtual void __Writeit(ostream &,string, int);
public:
	Gyration(){};
	Gyration(double a,Dvect & b, Dvect & c, Dvect & d): Radg(a),I(b),G(c), axis(d) {};
	Gyration(double a,Dvect & b, Dvect & c, Dvect & d, map<string,size_t> & hashtag): Radg(a),I(b),G(c), axis(d), str_hash(hashtag) {};
	Gyration(double a,Dvect & b, Dvect & c, Dvect & d, Dvect & xcm, map<string,size_t> & hashtag): Radg(a),I(b),G(c), axis(d), xcmx{new Dvect{xcm}}, str_hash(hashtag) {};
	virtual ~Gyration();
	double gRadg(){return Radg;}
	static json & gJson(){return myJson;}
	Dvect gI(){return I;}
	Dvect gG(){return G;}
	Dvect gaxis(){return axis;}
	Dvect gXcm(){return *xcmx;}
	void operator()(double a,Dvect & b, Dvect & c, Dvect & d){
		Radg=a;I=b;G=c;axis=d;
	}
	void operator()(double a,Dvect & b, Dvect & c, Dvect & d, Dvect & xcm){
		Radg=a;I=b;G=c;axis=d;
		if(!xcmx) xcmx=new Dvect{xcm};
		else *xcmx=xcm;
	}
	void operator()(double a,Dvect & b, Dvect & c, Dvect & d, double tt){
		Radg=a;I=b;G=c;axis=d;time_c=tt;
	}
	void operator()(double a,Dvect & b, Dvect & c, Dvect & d, Dvect & xcm, double tt){
		Radg=a;I=b;G=c;axis=d;time_c=tt;
		if(!xcmx) xcmx=new Dvect{xcm};
		else *xcmx=xcm;
	}
	Gyration & operator=(const Gyration & y){
		Radg=y.Radg;
		I=y.I;
		G=y.G;
		axis=y.axis;
		return *this;
	}
	static void setTime(double tt){time_c=tt;};
	Gyration operator/(int n){
		Gyration tmp;
		tmp.Radg=Radg/static_cast<double>(n);
		tmp.I=I/=static_cast<double>(n);
		tmp.G=G/static_cast<double>(n);
		tmp.axis=axis/static_cast<double>(n);
		return tmp;
	}
	Dvect Axis(){Dvect Ax=axis;for(int o=0;o<DIM;o++) Ax[o]=sqrt(Ax[o]);return Ax;}
	double Volume(){return (4.0/3.0)*M_PI*sqrt(axis[XX]*axis[YY]*axis[ZZ]);}
	void operator=(double y){
		Radg=y;
		I=y;G=y;axis=0.0;
	}
	Gyration & operator+=(Gyration & y){
		Radg+=y.Radg;
		I+=y.I;G+=y.G;axis+=y.axis;
		return *this;
	}
	template <typename G>
	friend ostream & operator<<(ostream & , Gyration<G> &);
	template <typename G>
	friend ostream & operator<<(ostream & , vector<Gyration<G>> &);
	template <typename G>
	friend ostream & operator<<(ostream & , Gyration<G> *);
	template <typename G>
	friend ostream & operator<<(ostream & , vector<Gyration<G> *> &);
};
#endif /* GYRATION_H_ */
