/*
 * GyrationJSON.h
 *
 *  Created on: Dec 21, 2017
 *      Author: marchi
 */

#ifndef MOLECULES_GYRATIONJSON_H_
#define MOLECULES_GYRATIONJSON_H_

#include "Gyration.h"
#include "json.hpp"
#include <set>
#include <string>

using std::set;
using std::string;
using nlohmann::json;
template <typename T>
class GyrationJSON: public Gyration<T> {
	static set<string> Time;
	void __Writeit(ostream &,string, int);
public:
	using Gyration<T>::Gyration;
	GyrationJSON(){}
	GyrationJSON(GyrationJSON<T> & y){
		this->Radg=y.Radg;
		this->I=y.I;
		this->G=y.G;
		this->axis=y.axis;
	}

	virtual ~GyrationJSON();
};

#endif /* MOLECULES_GYRATIONJSON_H_ */
