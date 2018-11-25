/*
 * DensMode.cpp
 *
 *  Created on: Nov 24, 2018
 *      Author: marchi
 */

#include "DensMode.h"
string DensMode::complement(string x){
	try{
		if(x == rSpace) return qSpace;
		if(x == qSpace) return rSpace;
		throw string("\nState "+ x +" is not a valid State. Only "+rSpace+
				" and "+qSpace +"are allowed\n");
	}catch(const string & s){
		cout << s << endl;
		exit(1);
	};
	return "";
}
void DensMode::operator()(string x){
	try{
		if(State.find(x) != State.end())
			State[x]=true;
		else throw string("\nState "+ x +" is not a valid State. Only "+rSpace+
				" and "+qSpace +"are allowed. Stop here.\n");
		State[complement(x)]=false;
	}catch(const string & s){
		cout << s << endl;
		exit(1);
	}

};
bool DensMode::operator[](string x){
	try{
		if(State.find(x) != State.end())
			return State[x];
		else throw string("\nState "+ x +" is not a valid State. Only "+rSpace+
				" and "+qSpace +"are allowed. Stop here.\n");
	}catch(const string & s){
		cout << s << endl;
		exit(1);
	}
}


