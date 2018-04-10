/*
 * ScatteringFactors.h
 *
 *  Created on: Jun 22, 2015
 *      Author: marchi
 */

#ifndef SRC_SCATTERINGFACTORS_H_
#define SRC_SCATTERINGFACTORS_H_
#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include "Finalize.h"

using namespace std;
/** \brief Handles the atomic scattering factors
 *         All factors as a function of q are store in a map of lambda functions
 *
 */
class ScatteringFactors {
public:
	struct opsfact{
		double a1{0},a2{0},a3{0},a4{0},b1{0},b2{0},b3{0},b4{0},c{0},d{1};
		string STR{""};
		opsfact(){};
		opsfact(const string);
		opsfact(const string, double);
		inline double operator()(double q1){
			double q=q1/(4.0*M_PI*10.0);
			double q2=q*q;
			return d*(a1*exp(-b1*q2)+a2*exp(-b2*q2)+a3*exp(-b3*q2)+a4*exp(-b4*q2)+c);
		}
		inline double operator()(){
			return a1;
		}

		// q is in A^{-1} and q1=4*Pi*q0  q1 and q0 are in nm!!
	};
	static bool hydrogens;
protected:
	static map<string,vector<vector<double> > > it1992;
	string Pick(const string s1) const{
		string s=s1;
		for(auto o=1;o<s1.size();o++)
			s[o]=tolower(s[o]);
		s[0]=toupper(s[0]);
//		if(s =="Na") s+="1+";
//		if(s =="Mg") s+="2+";
//		if(s =="Ca") s+="2+";
		double Fact=0.0;
		try{
			if(it1992.find(s) == it1992.end()) throw "Structure factor not found for "+ s;
		} catch(const string & ss){
			cout << ss << endl;
			Finale::Finalize::Final();
		}
		return s;
	}
	map<const string,opsfact> mySFact;
public:
	ScatteringFactors();

	const vector<vector<double> > & operator[](const string & s) const {
		return it1992.at(s);
	}
	double structFactDirect(const string & , double );
	double structFactDirect(const string & );
	const opsfact & operator()(const string S, double b=1.0){
		return this->StructFactor(S,b);
	}

	const opsfact & StructFactor(const string, double b=1.0);

	virtual ~ScatteringFactors();
};

#endif /* SRC_SCATTERINGFACTORS_H_ */
