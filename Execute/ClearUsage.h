/*
 * ClearUsage.h
 *
 *  Created on: Dec 20, 2017
 *      Author: marchi
 */

#ifndef SRC_CLEARUSAGE_H_
#define SRC_CLEARUSAGE_H_
#include <iterator>
#include <initializer_list>
#include <vector>
#include <map>
#include <string>
#include <iostream>

using std::vector;
using std::map;
using std::string;
using std::cout;
using std::endl;
class ClearUsage {
	vector<int> List;
	ClearUsage();
public:
	ClearUsage(std::initializer_list<int> list){for(auto op:list) List.push_back(op);};
	void ClrU(map<int,string> & Usage){
		for(auto op: List) Usage[op].clear();
	}
	virtual ~ClearUsage();
};

#endif /* SRC_CLEARUSAGE_H_ */
