/*
 * AtomsProp.h
 *
 *  Created on: Mar 1, 2019
 *      Author: marchi
 */

#ifndef PROPERTIES_ATOMSPROP_H_
#define PROPERTIES_ATOMSPROP_H_
#include "Atoms.h"
#include "properties.hpp"
#include "RhoHistogram.h"
#include <map>

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

public:
	AtomsProp(): Atoms<T>(){}
	AtomsProp(const int n): Atoms<T>(n){};
	AtomsProp(const AtomIndex & x):Atoms<T>(x){};
	void doProperty();
	void printProperty();
	void doTest(){cout << "Test is fine !" <<endl;};
};
template <typename T>
class AtomsProp<T,noprop>: public Atoms<T>{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;

public:
	AtomsProp(): Atoms<T>(){}
	AtomsProp(const int n): Atoms<T>(n){};
	AtomsProp(const AtomIndex & x):Atoms<T>(x){};
	void doProperty(){}
	void doTest(){cout << "Test is not fine !" <<endl;};
};

#endif /* PROPERTIES_ATOMSPROP_H_ */
