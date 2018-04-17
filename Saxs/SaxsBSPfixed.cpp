/*
 * SaxsBSPfixed.cpp
 *
 *  Created on: Feb 8, 2018
 *      Author: marchi
 */

#include "SaxsBSPfixed.h"

void SaxsBSPfixed::ComputeSAXS(RhoSaxs * Rho_ex,const MAtoms * y){
	AtomsD * x0{new MAtoms(*y)};

	if(!x0->CenterAtoms()) return;
/*
 * Pick padding. Padding is either average padding from box edges or custom for all systems except
 * for solute-only systems. In that case only zero padding is considered
 */
	Enums::Padding myPadding=pickPadding(y->getTypeNo());
	if(myPadding != zero) myPadding=Rho_ex->WhichPadding();
	void (RhoSaxs::*Padding)(RhoSaxs *,int,int,int){nullptr};
/*
 * Pick a pointer-to-member function to apply padding according to user's choice
 */
	switch(myPadding){
	case zero:
		cout << "zero"<<endl;
		Padding=&RhoSaxs::doPadding<zero>;
		break;
	case avgDensity:
		Padding=&RhoSaxs::doPadding<avgDensity>;
		break;
	case myDensity:
		Padding=&RhoSaxs::doPadding<myDensity>;
		break;
	case Periodic:
		Padding=&RhoSaxs::doPadding<Periodic>;
		break;
	default:
		break;
	}

	vector<vector<int> > g=y->getSaxsSolute();
	int tot{0};
	for(auto op: g)for(auto oq:op)tot++;
	Nsolute=tot;
	RhoSaxs & Rho_e=*Rho_ex;
	Nx=Rho_ex->getnnx();
	Ny=Rho_ex->getnny();
	Nz=Rho_ex->getnnz();

	if(!time){
		Rho_est=vector<RhoSaxs>(Sfacts.size(),RhoSaxs(Nx,Ny,Nz));
		Rho_fake=new RhoSaxs(Nx,Ny,Nz);
	}
	RhoSaxs & Rho_f=*Rho_fake;
	Mt=x0->getMt();
	Matrix oc=Mt.getOC();
	Matrix co=Mt.getCO();
	MCO+=co;
	MOC+=oc;
	Vol+=Mt.getVol();
	size_t nfx,nfy,nfz;
	nfx=(nx % 2 == 0)? nx/2: nx/2+1;
	nfy=(ny % 2 == 0)? ny/2: ny/2+1;
	nfz=(nz % 2 == 0)? nz/2: nz/2+1;
	Complex vt,vt2,v0;
	double Nt=static_cast<double> (Ntot);
	for(auto o=0;o<DIM;o++) {co[o][o]=SuperCell;oc[o][o]=1.0/SuperCell;}
	map<const string,ScatteringFactors::opsfact>::iterator it;
	this->__shift(x0);
	int Nsfact{0};
	map<string,double> MyPadding=Rho_ex->getPadding();
	for(it=Sfacts.begin();it!=Sfacts.end();it++){
		Rho_e=0.0;
		Rho_e.Density(order,x0,iSfacts[it->first],it->first);
		Rho_est[Nsfact]+=Rho_e;
		Rho_f+=Rho_e;
		cout << it->first<<endl;
		Nsfact++;
	}
	Rho_f.MakeAvg();
//	for (auto pad: Rho_ex->getPadding())
//		cout << pad.first<<" " << pad.second <<endl;
//	cout << Nx << " " << nx <<endl;
//	exit(1);
	time++;
	if(!(time%nskip)){
		int N_t=0;
		RhoSaxs Rho_alt(nx,ny,nz);
		array3<Complex> ro_k(nx,ny,nzp,align);
		array3<Complex> ro_ktot(nx,ny,nzp,align);
		array4<double> & ro_r=Rho_alt;
		crfft3d Backward3(nx,ro_k,ro_r[0]);
		rcfft3d Forward3(nx,ro_r[0],ro_k);

		int Nsfact{0};
		for(it=Sfacts.begin();it!=Sfacts.end();it++){
			Rho_est[Nsfact]/=static_cast<double> (nskip);
			(Rho_alt.*Padding)(&Rho_est[Nsfact],Nx,Ny,Nz);

			ro_k=Complex{0.0,0.0};
			Forward3.fft(ro_r[0],ro_k);
			auto ff=it->second;
			auto Na=iSfacts[it->first].size();
			N_t+=Na;
			Rho_est[Nsfact]=0.0;
			Nsfact++;
#pragma omp parallel for
			for(auto i=0;i<nx;i++){
				int ia=(i<nfx)?i : i-nx;
				for(auto j=0;j<ny;j++){
					int ja=(j<nfy)?j : j-ny;
					for(auto k=0;k<nzp;k++){
						int ka=(k<nfz)?k : k-nz;
						double mw1,mw2,mw3,mw;
						mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
						mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
						mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
						mw1=2.0*M_PI*mw1;
						mw2=2.0*M_PI*mw2;
						mw3=2.0*M_PI*mw3;
						mw=sqrt(mw1*mw1+mw2*mw2+mw3*mw3);
						double fq=ff(mw);
						ro_ktot[i][j][k]+=fq*ro_k[i][j][k];
					}
				}
			}
		}
		Modulus(ro_ktot,ro_ktot);
		Complex iD{1.0/N_t,0.0};
		ro_ktot*=iD;
		I_k+=ro_ktot;
		count++;
	}
	delete x0;
}

SaxsBSPfixed::~SaxsBSPfixed() {
	// TODO Auto-generated destructor stub
}

