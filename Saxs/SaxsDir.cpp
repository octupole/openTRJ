/*
 * SaxsDir.cpp
 *
 *  Created on: Mar 19, 2016
 *      Author: marchi
 */

#include "SaxsDir.h"

SaxsDir::SaxsDir() {
	// TODO Auto-generated constructor stub

}
void SaxsDir::ComputeSF(RhoSaxs * dummy,const MAtoms * y){
	const MAtoms & x=*y;
	qdfx=new SaxsHistogram(dq,qcut,unitsQ);
	SaxsHistogram & rdq=*qdfx;

	const double LARGE{200};
	const double HX{0.005};
	size_t nHist=LARGE/HX+1;
	vector<double> rdf(nHist);
	for(auto it0=Sfacts.begin();it0!=Sfacts.end();it0++){
		auto ff0=it0->second;
		vector<size_t> & ind0=iSfacts[it0->first];
		vector<Dvect> xa(ind0.size());
		double Na=ind0.size();
		for (auto o = 0; o < ind0.size(); o++) {
			size_t i=ind0[o];
			Dvect xx=Dvect{x.getX()[i]};
			xa[o][XX] = xx[XX];
			xa[o][YY] = xx[YY];
			xa[o][ZZ] = xx[ZZ];
		}
		for(auto it1=it0;it1!=Sfacts.end();it1++){
			auto ff1=it1->second;
			vector<size_t> & ind1=iSfacts[it1->first];
			double Nb=ind1.size();
			vector<Dvect> xb(ind1.size());
			for (auto o = 0; o < ind1.size(); o++) {
				size_t i=ind1[o];
				Dvect xx=Dvect{x.getX()[i]};
				xb[o][XX] = xx[XX];
				xb[o][YY] = xx[YY];
				xb[o][ZZ] = xx[ZZ];
			}
			double fact=it1==it0?1.0:2.0;
			for(auto & op: rdf) op=0.0;
			for(auto o=0;o<ind0.size();o++){
				Dvect xo=xa[o];
				for(auto p=0;p<ind1.size();p++){
					Dvect xx=xb[p]-xo;
					double rs=xx.Norm();
					if(rs > LARGE) continue;
					int h0=static_cast<int>(rs/HX);
					int h1=static_cast<int>(rint(rs/HX));
					if(h0 >= nHist || h1 >= nHist) exit(1);
					rdf[h0]+=fact*0.5;
					rdf[h1]+=fact*0.5;
				}

			}
			cout << " Doing interaction "<<it0->first + " "+ it1->first <<endl;
			for(auto o=0;o<rdq.HisX;o++){
				double qq=dq*o;

				double sum{0};
				for(auto w=0;w<rdf.size();w++){
					double r=w*HX;
					double qr=qq*r;
					sum+=qr==0?rdf[w]:rdf[w]*sin(qr)/qr;
				}

				double ff=ff0(qq)*ff1(qq);
				rdq[o]+=hist1D(ff*sum,0);
			}
		}
	}
	for(auto o=0;o<rdq.HisX;o++) rdq[o]=hist1D(rdq[o].hisX,1);
}

SaxsDir::~SaxsDir() {
	// TODO Auto-generated destructor stub
}

