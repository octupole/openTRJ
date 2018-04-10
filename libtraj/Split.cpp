/*
 * Split.cpp
 *
 *  Created on: Aug 8, 2017
 *      Author: marchi
 */

#include "../libtraj/Split.h"
const string toErase{"\t \n"};
vector<string> split(const string & s){
	std::stringstream iss(s);
	vector<string> tokens;
	copy(std::istream_iterator<string>(iss),
			std::istream_iterator<string>(),
			std::back_inserter<vector<string> >(tokens));
	for(unsigned int o=0;o<tokens.size();o++)
		tokens[o].erase(remove_if(tokens[o].begin(),tokens[o].end(),::isspace),tokens[o].end());
	return tokens;
}
void cleanString(string & s){
	std::size_t found = s.find_first_of(toErase);
	while(found != string::npos){
		s.erase(found,1);
		found=s.find_first_of(toErase,found+1);
	}
}
