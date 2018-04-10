/*
 * VectRotFlip.cpp
 *
 *  Created on: Feb 19, 2016
 *      Author: marchi
 */

#include "VecRotate.h"

VecRotate::VecRotate(vector<double> & y,int N) {
	try{
		if(y.size()%N) {
			stringstream ss0,ss1;
			ss0<<N;
			ss1<<y.size();
			throw string("\nSize of vector ("+ss1.str()+") must be a multiple of the rotation length ("+ss0.str()+").\n");
		}
	}catch(const string & s){cout << s<<endl;exit(1);}
	vT=y;
	Nx=N;
}
vector<double> VecRotate::vTget(int M){
	try{
		if(M+1 > Nx) throw string(" Improper use of vTget!");
	}catch(const string & s){cout << s <<endl;exit(1);}
	int Nbin=vT.size()/Nx;
	int nRot=vT.size()-Nbin*M;
	vector<double> vvT=vT;
	if(!nRot) return vvT;
	rotate(vvT.begin(), vvT.begin()+nRot,vvT.end());
	return vvT;
}
vector<double> VecFlipRotate::vTget(int M){
	try{
		if(M+1 > Nx) throw string(" Improper use of vTget!");
	}catch(const string & s){cout << s <<endl;exit(1);}
	int Nbin=vT.size()/Nx;
	M=Nx-M-1;

	int nRot=Nbin*M;
	vector<double> vvT=vT;
	std::reverse(vvT.begin(),vvT.end());
	if(!nRot) return vvT;
	rotate(vvT.begin(), vvT.begin()+nRot,vvT.end());
	return vvT;
};



