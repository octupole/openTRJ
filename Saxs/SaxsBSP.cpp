/*
 * SaxsBSP.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: marchi
 */

#include "SaxsBSP.h"

SaxsBSP::SaxsBSP():Saxs::Saxs() {}
SaxsBSP::SaxsBSP(int MyOrder, double Mydq, double Myqcut):Saxs::Saxs(MyOrder,Mydq,Myqcut) {
	BSpline::SetOrder(MyOrder);
}
void SaxsBSP::__shift(AtomsD * atm){
	AtomsD & x=*atm;
	double dx=0.5*(double)order/(double) Nx;
	double dy=0.5*(double)order/(double) Ny;
	double dz=0.5*(double)order/(double) Nz;
	MetricD Mt=atm->getMt();
	Matrix oc{Mt.getOC()};
	Matrix co{Mt.getCO()};
	vector<Dvect> xa=atm->getXA();
	for (auto i = 0; i < atm->getNR(); i++) {
		xa[i][XX] = oc[XX][XX] * x[i][XX] + oc[XX][YY] * x[i][YY] + oc[XX][ZZ] * x[i][ZZ]+dx;
		xa[i][YY] = oc[YY][XX] * x[i][XX] + oc[YY][YY] * x[i][YY] + oc[YY][ZZ] * x[i][ZZ]+dy;
		xa[i][ZZ] = oc[ZZ][XX] * x[i][XX] + oc[ZZ][YY] * x[i][YY] + oc[ZZ][ZZ] * x[i][ZZ]+dz;
		x[i][XX] =co[XX][XX]*xa[i][XX]+co[XX][YY]*xa[i][YY]+co[XX][ZZ]*xa[i][ZZ];
		x[i][YY] =co[YY][XX]*xa[i][XX]+co[YY][YY]*xa[i][YY]+co[YY][ZZ]*xa[i][ZZ];
		x[i][ZZ] =co[ZZ][XX]*xa[i][XX]+co[ZZ][YY]*xa[i][YY]+co[ZZ][ZZ]*xa[i][ZZ];
	}
}
void SaxsBSP::Setup(const vector<string> & at,bool mySans){
	Saxs::Setup(at,mySans);
	bsp_modx=new BSpmod(nx,ny,nz);
}
void SaxsBSP::Setup(const vector<int> & Lst,const vector<string> & at,bool mySans){
	Saxs::Setup(Lst,at,mySans);
	bsp_modx=new BSpmod(nx,ny,nz);
}
void SaxsBSP::Modulus(array3<Complex> & ro_k,array3<Complex> & ro_k1){
	size_t nzp=nz/2+1;
	double mw1,mw2,mw3,mw;
	size_t nfx,nfy,nfz;
	nfx=(nx % 2 == 0)? nx/2: nx/2+1;
	nfy=(ny % 2 == 0)? ny/2: ny/2+1;
	nfz=(nz % 2 == 0)? nz/2: nz/2+1;
	Matrix oc=Mt.getOC();
	for(auto i=0;i<nx;i++){
		double bsp_i=bsp_modx->ModuliX(i);
		for(auto j=0;j<ny;j++){
			double bsp_j=bsp_modx->ModuliY(j);
			for(auto k=0;k<nzp;k++){
				double bsp=bsp_i*bsp_j*bsp_modx->ModuliZ(k);
				ro_k[i][j][k]*=conj(ro_k1[i][j][k])*bsp;
			}
		}
	}

}


SaxsBSP::~SaxsBSP() {
	// TODO Auto-generated destructor stub
	delete bsp_modx;
}

