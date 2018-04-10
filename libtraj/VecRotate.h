/*
 * VectRotFlip.h
 *
 *  Created on: Feb 19, 2016
 *      Author: marchi
 */

#ifndef VECROTATE_H_
#define VECROTATE_H_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <functional>


using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::stringstream;
using std::rotate;

class VecRotate {
protected:
	vector<double> vT;
	int Nx{0};
public:
	VecRotate()=delete;
	VecRotate(vector<double>&,int);
	virtual vector<double> vTget(int);
	virtual ~VecRotate(){};
};

class VecFlipRotate: public VecRotate {
public:
	VecFlipRotate()=delete;
	VecFlipRotate(vector<double>& x,int y):VecRotate(x,y){};
	virtual vector<double> vTget(int);
	virtual ~VecFlipRotate(){};
};
#endif /* VECROTATE_H_ */
