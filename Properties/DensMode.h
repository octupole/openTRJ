/*
 * DensMode.h
 *
 *  Created on: Nov 24, 2018
 *      Author: marchi
 */

#ifndef EXECUTE_SAXS_DENSMODE_H_
#define EXECUTE_SAXS_DENSMODE_H_
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <ostream>
using std::map;
using std::cout;
using std::endl;
using std::string;

class DensMode {
	string rSpace{"R"};
	string qSpace{"Q"};
	std::map<string,bool> State{{rSpace,true},{qSpace,false}};
	string complement(string);
public:
	void operator()(string );
	bool operator[](string);
	friend std::ostream & operator<<(std::ostream & fout, DensMode & y){
		fout << "state R "<< y.State["R"]<< " state Q "<< y.State["Q"]<<endl;
		return fout;
	}
};

#endif /* EXECUTE_SAXS_DENSMODE_H_ */
