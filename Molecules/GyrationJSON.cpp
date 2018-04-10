/*
 * GyrationJSON.cpp
 *
 *  Created on: Dec 21, 2017
 *      Author: marchi
 */

#include "GyrationJSON.h"
template <typename T>
set<string> GyrationJSON<T>::Time;


template <typename T>
void GyrationJSON<T>::__Writeit(ostream & fout, string label, int o){
	static bool firstTime{true};
	string TimeC=std::to_string(this->time_c);
	if(!Time.count(TimeC)){
		Time.insert(TimeC);
		if(firstTime){
			firstTime=false;
			fout <<"\""<<TimeC<<"\": ";
		}else{
			fout <<"\""<<"gyro"<<"\": ";
			fout<< this->myJson;
			this->myJson.clear();
			fout <<"},";
			fout <<"\""<<TimeC<<"\": ";
		}

	}
	auto mySqrt=this->mySqrt;

	json & myType=this->myJson[label];
	json myClust;
	myClust["Rg"]=mySqrt(this->Radg);
	myClust["I"].push_back(mySqrt(this->I[XX]));
	myClust["I"].push_back(mySqrt(this->I[YY]));
	myClust["I"].push_back(mySqrt(this->I[ZZ]));
	myClust["ax"].push_back(mySqrt(this->axis[XX]));
	myClust["ax"].push_back(mySqrt(this->axis[YY]));
	myClust["ax"].push_back(mySqrt(this->axis[ZZ]));
	myClust["xcm"].push_back((*this->xcmx)[XX]);
	myClust["xcm"].push_back((*this->xcmx)[YY]);
	myClust["xcm"].push_back((*this->xcmx)[ZZ]);
	for(auto tags: this->str_hash){
		myClust["hsh"][tags.first]=tags.second;
	}
	myType.push_back(myClust);
}

template <typename T>
GyrationJSON<T>::~GyrationJSON() {
	// TODO Auto-generated destructor stub
}

template class GyrationJSON<float>;
template class GyrationJSON<double>;
