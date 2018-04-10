/*
 * RhoSaxsBSP.cpp
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#include "RhoSaxsBSP.h"


void RhoSaxsBSP::__Density(const int pq,const AtomsD * y, vector<size_t> & ind, string myType
		,vector<double> & wei){
	Type=myType;
	AtomsD  x{*y};
	int order=BSpline::Order();
	int nx0=static_cast<int>(nnx);
	int ny0=static_cast<int>(nny);
	int nz0=static_cast<int>(nnz);
	size_t natoms=ind.size();
	vector<Dvect> xa=vector<Dvect>(natoms,0.0);
	Real x1,y1,z1,r1,s1,t1,gx,gy,gz;
	int mx,my,mz;
	int ox,oy,oz,tx,ty,tz;
	MetricD Mt=x.getMt();

	setMetric(Mt.getCO());
	double DV=RhoSaxs::getDV();
	for (auto o = 0; o < ind.size(); o++) {
		size_t i=ind[o];
		xa[o][XX] = oc[XX][XX] * x[i][XX] + oc[XX][YY] * x[i][YY] + oc[XX][ZZ] * x[i][ZZ];
		xa[o][YY] = oc[YY][XX] * x[i][XX] + oc[YY][YY] * x[i][YY] + oc[YY][ZZ] * x[i][ZZ];
		xa[o][ZZ] = oc[ZZ][XX] * x[i][XX] + oc[ZZ][YY] * x[i][YY] + oc[ZZ][ZZ] * x[i][ZZ];
	}

	BSpline MySplineX;
	BSpline MySplineY;
	BSpline MySplineZ;

	for (auto s = 0; s < ind.size(); s++) {
		auto chg=wei[s];
		x1 = xa[s][XX];
		y1 = xa[s][YY];
		z1 = xa[s][ZZ];
		r1 = static_cast<Real> (nx0 * (x1 - rint(x1 - 0.5)));
		s1 = static_cast<Real> (ny0 * (y1 - rint(y1 - 0.5)));
		t1 = static_cast<Real> (nz0 * (z1 - rint(z1 - 0.5)));
		mx = static_cast<int>(r1);
		my = static_cast<int>(s1);
		mz = static_cast<int>(t1);
		gx = r1 - static_cast<double>(mx);
		gy = s1 - static_cast<double>(my);
		gz = t1 - static_cast<double>(mz);
		spline splX=MySplineX(gx);
		spline splY=MySplineY(gy);
		spline splZ=MySplineZ(gz);

		int i0=mx-order;
		for(auto  o=0;o<order;o++){
			double fact_o=chg*splX.x[o];
			int i=i0+(nx0-((i0>=0)?nx0:-nx0))/2;
			int j0=my-order;
			for(auto  p=0;p<order;p++){
				double fact_p=fact_o*splY.x[p];
				int j=j0+(ny0-((j0>=0)?ny0:-ny0))/2;
				int k0=mz-order;
				for(auto q=0;q<order;q++){
					double fact_q=fact_p*splZ.x[q];
					int k=k0+(nz0-((k0>=0)?nz0:-nz0))/2;
					(*this)[0][i][j][k]+=fact_q;
					k0++;
				}
				j0++;
			}
			i0++;
		}
	}


}
RhoSaxsBSP::~RhoSaxsBSP() {
	// TODO Auto-generated destructor stub
}

