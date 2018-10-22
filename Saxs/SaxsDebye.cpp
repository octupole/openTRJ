/*
 * SaxsDebey.cpp
 *
 *  Created on: Mar 19, 2016
 *      Author: marchi
 */

#include "SaxsDebye.h"

SaxsDebye::SaxsDebye() {
	// TODO Auto-generated constructor stub

}

void SaxsDebye::ComputeSAXS(RhoSaxs * dummy,const MAtoms * y){

	const MAtoms & x=*y;
	MetricD Mt=x.getMt();
	Matrix OC{Mt.getOC()};
	Matrix CO{Mt.getCO()};
	mycut=0.5*(CO[XX][XX]+CO[YY][YY]+CO[ZZ][ZZ])/3.0;
	qdfx=new SaxsHistogram(dq,qcut,unitsQ);
	SaxsHistogram & rdq=*qdfx;

	for(auto q=0;q<rdq.HisX;q++){
		double qq=dq*q;
		double sum=0.0;
		cout << " Doing q = " << q*dq <<endl;
		for(auto it0=Sfacts.begin();it0!=Sfacts.end();it0++){
			auto ff0=it0->second;
			vector<size_t> & ind0=iSfacts[it0->first];
			vector<Dvect> xa(ind0.size());
			for (auto o = 0; o < ind0.size(); o++) {
				size_t i=ind0[o];
				Dvect xx=Dvect{x.getXA()[i]};
				xa[o][XX] = xx[XX];
				xa[o][YY] = xx[YY];
				xa[o][ZZ] = xx[ZZ];
			}
			for(auto it1=it0;it1!=Sfacts.end();it1++){
				auto ff1=it1->second;
				auto ff=ff0(qq)*ff1(qq);
				vector<size_t> & ind1=iSfacts[it1->first];
				vector<Dvect> xb(ind1.size());
				for (auto o = 0; o < ind1.size(); o++) {
					size_t i=ind1[o];
					Dvect xx=Dvect{x.getXA()[i]};
					xb[o][XX] = xx[XX];
					xb[o][YY] = xx[YY];
					xb[o][ZZ] = xx[ZZ];
				}
				double fact=it1==it0?1.0:2.0;
				for(auto o=0;o<ind0.size();o++){
					Dvect xo=xa[o];
					for(auto p=0;p<ind1.size();p++){
						Dvect xx=xb[p]-xo;
						for(int q=0;q<DIM;q++) xx[q]=xx[q]-rint(xx[q]);
						Dvect xc=CO*xx;
						double x=xc.Norm();
						double val=qq*x==0.0?1.0:sin(qq*x)/(qq*x);
						sum+=fact*ff*val;
					}

				}
			}
		}
		rdq[q]+=hist1D{sum,1};
	}
}

SaxsDebye::~SaxsDebye() {
	// TODO Auto-generated destructor stub
}

