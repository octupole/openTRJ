/*
 * AtomsProp.h
 *
 *  Created on: Mar 1, 2019
 *      Author: marchi
 */

#ifndef PROPERTIES_ATOMSPROP_H_
#define PROPERTIES_ATOMSPROP_H_
#include "Atoms.h"
#include <cmath>
#include <map>
#include <tuple>

using std::map;

template <typename T, myOptions OPT>
class AtomsProp;

template <typename T>
class AtomsProp<T,radial>: public Atoms<T>{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	map<string,Properties::RhoHistogram *> histograms;
	map<string,double> masses;
	static string headerXVG;
	T rcut{0};

public:
	AtomsProp(): Atoms<T>(){}
	AtomsProp(const int n): Atoms<T>(n){};
	AtomsProp(const AtomIndex & x):Atoms<T>(x){};
	void * doProperty();
	myOptions getProperty(){return myOptions::radial;}
	void Reduce(Parallel::NewMPI *);
};
template <typename T>
class AtomsProp<T,gyro>: public Atoms<T>{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	std::tuple<vector<Gyration<T>*>, Percolation<T> * > out;

public:
	AtomsProp(): Atoms<T>(){}
	AtomsProp(const int n): Atoms<T>(n){};
	AtomsProp(const AtomIndex & x):Atoms<T>(x){};
	void * doProperty();
	myOptions getProperty(){return myOptions::gyro;}
	void Reduce(Parallel::NewMPI *){};

};
template <typename T>
class AtomsProp<T,gyroJ>: public Atoms<T>{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	std::tuple<vector<Gyration<T>*>, Percolation<T> *> out;

public:
	AtomsProp(): Atoms<T>(){}
	AtomsProp(const int n): Atoms<T>(n){};
	AtomsProp(const AtomIndex & x):Atoms<T>(x){};
	void * doProperty();
	myOptions getProperty(){return myOptions::gyroJ;}
	void Reduce(Parallel::NewMPI *){};

};

#endif /* PROPERTIES_ATOMSPROP_H_ */
