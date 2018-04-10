/*
 * RhoSaxs.cpp
 *
 *  Created on: Aug 19, 2015
 *      Author: marchi
 */

#include <RhoSaxs.h>
bool RhoSaxs::firstTime=true;
void RhoSaxs::Density(const int pq, const AtomsD * y, vector<size_t> & ind, string myType){
	vector<double> wei(ind.size(),1.0);
	__Density(pq,y,ind,myType,wei);
}
void RhoSaxs::Density(const int pq, const AtomsD * y, vector<size_t> & ind, string myType
		, vector<double> & wei){
	try{
		if(wei.size()!= ind.size()) throw string("Weights for density have wrong dimensions.");
	}catch(const string & s){cout << s<<endl;exit(1);}
	__Density(pq,y,ind,myType,wei);
}

void RhoSaxs::__Density(const int pq, const AtomsD * y, vector<size_t> & ind, string myType
		, vector<double> & wei){
	Type=myType;
	AtomsD  x{*y};
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
	for (unsigned int i = 0; i < ind.size(); i++) {
		Real chg[1];
		chg[0]=static_cast<Real>(wei[i]);
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
		TriLinear(gx, gy, gz, chg);
		for (int n = 0; n < VERTEX; n++) {
			tx = static_cast<int>(cube[n][XX]);
			ty = static_cast<int>(cube[n][YY]);
			tz = static_cast<int>(cube[n][ZZ]);
			ox = (nx0 > mx + tx) ? mx + tx : mx + tx - nx0;
			oy = (ny0 > my + ty) ? my + ty : my + ty - ny0;
			oz = (nz0 > mz + tz) ? mz + tz : mz + tz - nz0;
			(*this)[0][ox][oy][oz] += static_cast<double>(qq[n][0]);
		}
	}
}
auto BorderBins=[](size_t nn, double shell)->size_t{
	return static_cast<size_t> (shell*nn/2);
};
void RhoSaxs::zeroPadding(RhoSaxs * y, int Nx, int Ny, int Nz){
	if(firstTime){
		firstTime=false;
		cout << "\nDoing FFT zero padding.\n" <<endl;
	}
	int ifx=nnx/2-Nx/2;
	int ify=nny/2-Ny/2;
	int ifz=nnz/2-Nz/2;
	for(auto o=0;o<Nx;o++)for(auto p=0;p<Ny;p++)for(auto q=0;q<Nz;q++)
		(*this)[0][o+ifx][p+ify][q+ifz]=(*y)[0][o][p][q];
	return;

}
void RhoSaxs::PBCPadding(RhoSaxs * y, int Nx, int Ny, int Nz){
	if(firstTime){
		firstTime=false;
		cout << "\nDoing FFT periodic padding.\n" <<endl;
	}
	for(auto o=0;o<nnx;o++)for(auto p=0;p<nny;p++)for(auto q=0;q<nnz;q++){
		auto _o=o%Nx;auto _p=p%Ny;auto _q=q%Nz;
		(*this)[0][o][p][q]=(*y)[0][_o][_p][_q];
	}


}
void RhoSaxs::avgPadding(RhoSaxs * y, int Nx, int Ny, int Nz){

	double N1=(double) nnx/(double) Nx;
	double N2=(double) nny/(double) Ny;
	double N3=(double) nnz/(double) Nz;
	if(firstTime){
		firstTime=false;
		cout << "\nDoing FFT padding with average border density\n" <<endl;
	}
// get border bins
	int Mx=BorderBins(Nx,SHELL);
	int My=BorderBins(Ny,SHELL);
	int Mz=BorderBins(Nz,SHELL);
	try{
		if(!Mx || !My || !Mz) throw string("The number of border bins in some direction is zero.\n ")+
				string("This happens if -lowq is large and the grid dimensions are small.");
	}catch(const string & s){cout <<s<<endl;Finale::Finalize::Final();}
	double summ0{0.0},summ1{0.0}, summ2{0};
	for(auto o=0;o<Nx;o++)for(auto p=0;p<Ny;p++)for(auto q=0;q<Nz;q++)
		summ0+=(*y)[0][o][p][q];
	for(auto o=Mx;o<Nx-Mx;o++)for(auto p=My;p<Ny-My;p++)for(auto q=Mz;q<Nz-Mz;q++)
		summ1+=(*y)[0][o][p][q];

	double myDens=(summ0-summ1)/(Nx*Ny*Nz-(Nx-2*Mx)*(Ny-2*My)*(Nz-2*Mz));
	(*this)=myDens;
//	cout << y->Type << " " <<myDens <<endl;

	for(auto o=0;o<Nx;o++)for(auto p=0;p<Ny;p++)for(auto q=0;q<Nz;q++)
		(*this)[0][o][p][q]=(*y)[0][o][p][q];

	summ0=-myDens*(N1*N2*N3-1.0)/(N1*N2*N3);
	for(auto o=0;o<nnx;o++)for(auto p=0;p<nny;p++)for(auto q=0;q<nnz;q++)
		(*this)[0][o][p][q]+=summ0;
}
void RhoSaxs::myPadding(RhoSaxs * y, int Nx, int Ny, int Nz){
	double N1=(double) nnx/(double) Nx;
	double N2=(double) nny/(double) Ny;
	double N3=(double) nnz/(double) Nz;
	if(firstTime){
		firstTime=false;
		cout << "\nDoing FFT padding with custom density\n" <<endl;
	}
	double vol=y->getVol();
	double No=vol*y->MyPadding[y->Type];
	double myDens=No/(Nx*Ny*Nz);

	(*this)=myDens;
	cout << y->Type << " " <<  y->MyPadding[y->Type]<<endl;

	double summ0{0.0};
	for(auto o=0;o<Nx;o++)for(auto p=0;p<Ny;p++)for(auto q=0;q<Nz;q++)
		(*this)[0][o][p][q]=(*y)[0][o][p][q];

	summ0=-myDens*(N1*N2*N3-1.0)/(N1*N2*N3);
	for(auto o=0;o<nnx;o++)for(auto p=0;p<nny;p++)for(auto q=0;q<nnz;q++)
		(*this)[0][o][p][q]+=summ0;
}
void RhoSaxs::selectPadding(Padding y){
	exPadding=y;
}

void RhoSaxs::setPadding(const map<string,double> & y){
	if(y.empty()) return;
	MyPadding=y;
}
RhoSaxs::~RhoSaxs() {
	// TODO Auto-generated destructor stub
}
void RhoSaxs::MakeAvg(){
	double avg{0};
	double * it=this[0][0][0][0];
	RhoSaxs RhoN(*this);
	RhoN=0.0;
	cout << RhoN.getnnx()<<endl;
	size_t nNZero{0};
	for(auto m=0u;m<nnx*nny*nnz;m++){
		double val=*(it++);
		if(val) {
			avg+=val;
			nNZero++;
		}
	}
	auto dens=avg/static_cast<double>(nNZero);
	it=this[0][0][0][0];
	for(auto m=0u;m<nnx*nny*nnz;m++){
		double val=*(it++);
		if(val) *it=dens;
	}
	cout << dens << endl;exit(1);
//	for(auto mx=0u;mx<nnx;mx++){
//		for(auto my=0u;my<nny;my++){
//			for(auto mz=0u;mz<nnz;mz++){
//				if( (*this)[0][mx][my][mz] ){
//					avg+=(*this)[0][mx][my][mz];
//				}
//			}
//		}

//	}
}
