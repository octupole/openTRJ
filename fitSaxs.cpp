/*
 ============================================================================
 Name        : FitSaxs.cpp
 Author      : 
 Version     :
Copyright :  Massimo Marchi Massimo.Marchi(at)cea.fr

This software is a computer program whose purpose is to fit experimental and
computed SAXS spectra.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 Description :
 ============================================================================

 */

#include <iostream>
#include <sstream>
#include <algorithm>
#include "TrjRead.h"
#include "stdafx.h"
#include "dataanalysis.h"
#include "Spline1DInterpolant.h"

using namespace alglib;

using namespace std;

auto ReadL=[](string str)->vector<double>{
	stringstream ss{str};
	double x,y,s;
	ss>>x>>y>>s;
	if(s == y) s=1.0;
	return vector<double>{x,y,s};
};
auto myComp=[](vector<double> x,vector<double> y) {return x[0]<y[0];};

int main(int argc, char ** argv)
{
	trj::TrjRead MyIn(argc,argv);
	MyIn.Input();
	ofstream & fout=*MyIn.gFoutx();
	ifstream & fin1=*MyIn.gFin1();
	ifstream & fin2=*MyIn.gFin2();
	vector<vector<double> > mp1,mp2;

	for(string str;getline(fin1,str);){
		if(str.find('#') != string::npos) continue;
		auto v=ReadL(str);
		mp1.push_back(v);
	}
	Spline1D::Spline1DInterpolant func1(mp1);
	for(string str;getline(fin2,str);){
		if(str.find('#') != string::npos) continue;
		auto v=ReadL(str);
		mp2.push_back(v);
	}
	std::sort(mp2.begin(),mp2.end(),myComp);
	real_2d_array xy;
	real_1d_array S;
	xy.setlength(mp2.size(),2);
	S.setlength(mp2.size());
	int o{0};
	for(auto op: mp2){
		xy[o][0]=op[1];
		xy[o][1]=func1(op[0]);
		S[o]=xy[o][1];
		o++;
	}
	ae_int_t info;
    ae_int_t nvars;
    linearmodel model;
    lrreport rep;
    real_1d_array c;
    lrbuilds(xy,S,mp2.size(), 1, info, model, rep);

    lrunpack(model, c, nvars);
    printf("%s\n", c.tostring(4).c_str());

    for(auto op: mp1){
    	double y=(op[1]-c[1])/c[0];
    	fout << op[0] <<" " << y <<endl;
    }
	return 0;
}
