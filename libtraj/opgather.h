/*
 * opgather.h
 *
 *  Created on: Jun 24, 2015
 *      Author: marchi
 */

#ifndef SRC_OPGATHER_H_
#define SRC_OPGATHER_H_
#include <string>
#include <vector>
using namespace std;

template<typename T>
class opgather {
	T ToCompare;
public:
	opgather(T y): ToCompare(y){};
	vector<size_t> operator()(const vector<T> & z){
		vector<size_t> data;
		size_t nat=0;
		for(auto it=z.begin();it != z.end();it++){
			if(*it == ToCompare) data.push_back(nat);
			nat++;
		}
		return data;
	}
	virtual ~opgather(){};
};

#endif /* SRC_OPGATHER_H_ */

