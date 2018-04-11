/*
 * TrjRead.h
 *
 *  Created on: Dec 23, 2015
 *      Author: marchi
 */


#ifndef SRC_TRJREAD_H_
#define SRC_TRJREAD_H_

#include "trjInput.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iterator>


using namespace std;

namespace trj {

//!
//! Create pointers to TrjRead class member streams
//!
template <typename T>
class Streams{
	T * myStream{nullptr};
public:
	Streams(){};
	Streams(T * y):myStream{&y}{};
	 Streams & operator=(T * y ){
		myStream=y;
		return *this;
	}
	T * operator()(){
		return myStream;
	}

};
//!
//! Create pointers to TrjRead class member values
//!
template <typename T>
class Values{
	T * const myValue;
public:
	Values(T & y): myValue{&y}{};
	T & operator()(){return *myValue;}
};

//! This class is used to create the input information from standard input
/*!
 *
 */
class TrjRead: public trjInput {

	string filein1,filein2,fileout="Fit.out";
	ifstream ftest;
	ofstream * foutx{nullptr};
	ifstream * fin1{nullptr};
	ifstream * fin2{nullptr};
public:
	Values<string> gfilein1{filein1};
	Values<string> gfilein2{filein2};
	Values<string> gfileout{fileout};
 	Streams<ofstream> gFoutx;
 	Streams<ifstream> gFin1;
	Streams<ifstream> gFin2;

	TrjRead(int nv,char ** v);//!< Constructor from input line
	void Input();
	virtual ~TrjRead();
};

} /* namespace trj */

#endif /* SRC_TRJREAD_H_ */
