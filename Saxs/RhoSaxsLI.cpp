/*
 * RhoSaxsLI.cpp
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#include "RhoSaxsLI.h"

RhoSaxsLI::RhoSaxsLI() {
	// TODO Auto-generated constructor stub

}
void RhoSaxsLI::__Density(const int order,const AtomsD * y, vector<size_t> & ind,string myType
		,vector<double> & wei){
	Type=myType;
	AtomsD x{*y};
	LagrangeInterpolation LI{order};
	auto  MyPoly=LI.Poly();
	vector<int> tvect;
	for(auto pv: MyPoly)
		tvect.push_back(pv.first);


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
	auto Idx=[&nx0](const int mx){return mx<0?mx+nx0:(mx<nx0?mx:mx-nx0);};
	auto Idy=[&ny0](const int my){return my<0?my+ny0:(my<ny0?my:my-ny0);};
	auto Idz=[&nz0](const int mz){return mz<0?mz+nz0:(mz<nz0?mz:mz-nz0);};

	for (unsigned int i = 0; i < ind.size(); i++) {
		auto chg=wei[i];
		x1 = xa[i][XX];
		y1 = xa[i][YY];
		z1 = xa[i][ZZ];
		r1 = static_cast<Real> (nx0 * (x1 - rint(x1 - 0.5)));
		s1 = static_cast<Real> (ny0 * (y1 - rint(y1 - 0.5)));
		t1 = static_cast<Real> (nz0 * (z1 - rint(z1 - 0.5)));
		mx = static_cast<int>(r1);
		my = static_cast<int>(s1);
		mz = static_cast<int>(t1);
		gx = r1 - (Real) mx;
		gy = s1 - (Real) my;
		gz = t1 - (Real) mz;
		for(auto o=0;o<tvect.size();o++){
			int tx=tvect[o]*chg;
			double ff_o=MyPoly[tx](gx);
			auto ox = Idx(mx + tx);

			for(auto p=0;p<tvect.size();p++){
				int ty=tvect[p];
				double ff_p=MyPoly[ty](gy)*ff_o;
				auto oy = Idy(my + ty);

				for(auto q=0;q<tvect.size();q++){
					int tz=tvect[q];
					double ff_q=MyPoly[tz](gz)*ff_p;
					auto oz = Idz(mz + tz);
					(*this)[0][ox][oy][oz] += ff_q;
				}
			}
		}

	}

}

RhoSaxsLI::~RhoSaxsLI() {
	// TODO Auto-generated destructor stub
}

