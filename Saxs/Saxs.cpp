/*
 * Saxs.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: marchi
 */

#include "Saxs.h"

const bool myDEBUG{true};
const double Avog{6.022140857e23};
Saxs::Saxs(){

}
void Saxs::setSuperCell(double y){
	SuperCell=y;
}
Saxs & Saxs::operator-=(Saxs & y){
	*qdfx=__Qhistogram();
	Spline1D::Spline1DInterpolant spline0(qdfx,dq,unitsQ);
	if(y.nx != 0 && y.ny != 0 && y.nz !=0){
		*y.qdfx=y.__Qhistogram();
		double vol0{1.0},vol1{1.0};
		for(auto o=0;o<DIM;o++){
			vol0*=MCO[o][o];
			vol1*=y.MCO[o][o];
		}
		double alpha=MassSolute>0?1.0-(MassSolute*1.0e24/(vol0*Avog))*0.743/1000.0:1.0;
		cout << "\n ---- Using subtraction factor " << alpha << " -----" <<endl;
		double Factor=vol1!=0.0?(vol0/vol1)*alpha:0.0;
		for(auto o=0;o<y.qdfx->HisX;o++){
			double value=(*y.qdfx)[o].Ratio()*Factor;
			(*y.qdfx)[o]=hist1D{value,1};
		}
		Spline1D::Spline1DInterpolant spline1(y.qdfx,y.dq,y.unitsQ);
		spline0-=spline1;
	}
	qdfx->clear();
	(*qdfx)(dq_orig,qcut);
	for(auto o=0;o<qdfx->HisX;o++){
		double diff=spline0(o*dq_orig);
		(*qdfx)[o]=hist1D{diff,1};
	}
	return *this;
}
void Saxs::__copy(const Saxs & y ){
	dq=y.dq; qcut=y.qcut; order=y.order;
	noSplineOut=y.noSplineOut;
	Allocate(y.nx,y.ny,y.nz);
	I_k=y.I_k;
	I_r=y.I_r;
	SuperCell=y.SuperCell;
	SuperCell0=y.SuperCell0;

	Nsolute=y.Nsolute;
	Ntot=y.Ntot;
	bAvg=y.bAvg;
	MCO=y.MCO;
	MOC=y.MOC;
	if(noSplineOut)
		qdfx=new SaxsHistogram;
	else
		qdfx=new SaxsHistogramSpline;

	*qdfx=*y.qdfx;
	count=y.count;
}
Saxs & Saxs::operator =(const Saxs & y){
	__copy(y);
	return *this;
}
Saxs::Saxs(const Saxs & y){
	__copy(y);
}
void Saxs::Allocate(size_t mx,size_t my,size_t mz){
	nx=mx;
	ny=my;
	nz=mz;
	try{
		if(nx == 0 || ny == 0 || nz == 0) throw string("Must define grid dimension first. Abort.");
	} catch(const string & s){
		cout << mx <<endl;
		cout << s <<endl;
		Finale::Finalize::Final();
	}
	nzp=nz/2+1;
	I_k.Allocate(nx,ny,nzp,align);
	I_r.Allocate(nx,ny,nz,align);
}
void Saxs::Reduce(Parallel::NewMPI * y){
	if(!y->Get_Size()) return;
	int ncount=this->count;
	y->ReduceSum(&ncount,1);
	this->count=ncount;
	auto dim=nx*ny*nzp;
	auto ip=&I_k[0][0][0];
	y->ReduceSum(ip,dim);

	dim=nx*ny*nz;
	auto ipr=&I_r[0][0][0];
	y->ReduceSum(ipr,dim);
	dim=DIM*DIM;
	auto ip2=&MCO[0][0];
	y->ReduceSum(ip2,dim);
	ip2=&MOC[0][0];
	y->ReduceSum(ip2,dim);
}

void Saxs::Averages(){
	bAvg=true;
	I_k/=Complex{static_cast<double>(count),0.0};
	I_r/=static_cast<double>(count);
	MCO/=static_cast<double>(count);
	MOC/=static_cast<double>(count);
}
SaxsHistogram Saxs::__Qhistogram(){
	Matrix oc{0.0};
	for(auto o=0;o<DIM;o++) {oc[o][o]=1.0/SuperCell;}
	size_t nfx{(nx % 2 == 0)? nx/2: nx/2+1},nfy{(ny % 2 == 0)? ny/2: ny/2+1}
		,nfz{(nz % 2 == 0)? nz/2: nz/2+1};

	double mw1,mw2,mw3,mw;
	Dvect fx{(double)nfx-1,(double)nfy-1,(double)nfz-1};
	double m_qcut{qcut},m_dq{dq};
	vector<double> mydq0={2.0*M_PI*oc[XX][XX],2.0*M_PI*oc[YY][YY],2.0*M_PI*oc[ZZ][ZZ],this->dq};
	vector<double> mycut0={2.0*M_PI*oc[XX][XX]*fx[XX],2.0*M_PI*oc[YY][YY]*fx[YY],2.0*M_PI*oc[ZZ][ZZ]*fx[ZZ],this->qcut};
	double dq=1.1*(*std::max_element(mydq0.begin(),mydq0.end()));
	double qcut=*std::min_element(mycut0.begin(),mycut0.end());
	try{
		int ntry{0};
		bool notok=false;
		stringstream ss;
		string msg("");
		if(dq != m_dq && noSplineOut){
			ss<< "from " <<m_dq*unitsQ << " to "<<dq*unitsQ ;
			msg+="    Histogram was constructed with a larger bin " + ss.str();
			notok=true;
			ntry++;
		}
		if(qcut != m_qcut) {
			ss.str(string());
			ss<< "from " <<m_qcut*unitsQ << " to "<<qcut*unitsQ ;
			if(ntry) msg+="\n";
			msg+="    Cutoff was too large. Changed to maximum allowed " + ss.str();
			notok=true;
			this->qcut=qcut;
			ntry++;
		}
		if(notok) throw msg;
	}catch(const string & s){
		cout << "\n***"+string(80,'-')+"***\n";
		cout << s <<endl;
		cout << "***"+string(80,'-')+"***\n";
	}
	SaxsHistogram qdf{dq,qcut,unitsQ};
	qdf.setDq(m_dq);
	for(auto i=0;i<nx;i++){
		int ia=(i<nfx)?i : i-nx;
		size_t ib=i==0?0:nx-i;
		for(auto j=0;j<ny;j++){
			int ja=(j<nfy)?j : j-ny;
			size_t jb=j==0?0:ny-j;
			for(auto k=0;k<nzp;k++){
				int ka=(k<nfz)?k : k-nz;
				mw1=oc[XX][XX]*ia+oc[XX][YY]*ja+oc[XX][ZZ]*ka;
				mw2=oc[YY][XX]*ia+oc[YY][YY]*ja+oc[YY][ZZ]*ka;
				mw3=oc[ZZ][XX]*ia+oc[ZZ][YY]*ja+oc[ZZ][ZZ]*ka;
				mw1=2.0*M_PI*mw1;
				mw2=2.0*M_PI*mw2;
				mw3=2.0*M_PI*mw3;
				mw=sqrt(mw1*mw1+mw2*mw2+mw3*mw3);
				if(mw<qcut){
					Complex v0=I_k[i][j][k];
					if(k != 0 && k != nzp-1){
						Complex vt1=I_k[ib][jb][k];
						v0=0.5*(v0+vt1);
					}
					int h0=static_cast<int>(mw/dq);
					int h1=h0+1;
					try{
						if(h0 >= qdf.HisX || h1 >= qdf.HisX) throw string("Something wrong here! Value outside histogram dimensions.");
					} catch(const string & s){cout << s <<endl;exit(1);}
					qdf[h0]+=hist1D(v0.real(),1);
					if(h0 != 0)
						qdf[h1]+=hist1D(v0.real(),1);


				}
			}
		}
	}
	for(auto h=0;h<qdf.HisX;h++) qdf[h]=hist1D(qdf[h].Ratio()*Ntot,1);
	return qdf;
}
void Saxs::GofR(){
	double mw1,mw2,mw3,mw;
	array3<Complex> ro_ktot(nx,ny,nzp,align);
	array3<double> ro_r(nx,ny,nz);
	crfft3d Backward3(nx,ro_ktot,ro_r);
	Matrix co=MCO;
	co[XX][XX]=SuperCell;
	co[YY][YY]=SuperCell;
	co[ZZ][ZZ]=SuperCell;
	size_t nfx{(nx % 2 == 0)? nx/2: nx/2+1},nfy{(ny % 2 == 0)? ny/2: ny/2+1}
		,nfz{(nz % 2 == 0)? nz/2: nz/2+1};
	for(auto i=0;i<nx;i++)
		for(auto j=0;j<ny;j++)
			for(auto k=0;k<nzp;k++)
				ro_ktot[i][j][k]=I_k[i][j][k];
	Backward3.fftNormalized(ro_ktot,ro_r);
	double dx{0.1},rcut{SuperCell};
	SaxsHistogram gofr{dx,rcut,unitsR};
	for(auto i=0;i<nx;i++){
		int ia=(i<nfx)?i : i-nx;
		double xa=static_cast<double>(ia)/static_cast<double>(nx);
		for(auto j=0;j<ny;j++){
			int ja=(j<nfy)?j : j-ny;
			double ya=static_cast<double>(ja)/static_cast<double>(ny);
			for(auto k=0;k<nzp;k++){
				int ka=(k<nfz)?k : k-nz;
				double za=static_cast<double>(ka)/static_cast<double>(nz);
				mw1=co[XX][XX]*xa+co[XX][YY]*ya+co[XX][ZZ]*za;
				mw2=co[YY][XX]*xa+co[YY][YY]*ya+co[YY][ZZ]*za;
				mw3=co[ZZ][XX]*xa+co[ZZ][YY]*ya+co[ZZ][ZZ]*za;
				mw=sqrt(mw1*mw1+mw2*mw2+mw3*mw3);
				if(mw<rcut){
					double v0=ro_r[i][j][k];
					int h0=static_cast<int>(mw/dx);
					try{
						if(h0 >= gofr.HisX) throw string("Something wrong here! Value outside histogram dimensions.");
					} catch(const string & s){cout << s <<endl;exit(1);}
					gofr[h0]+=hist1D(v0,1);
				}
			}
		}
	}
	for(auto h=0;h<gofr.HisX;h++) gofr[h]=hist1D(gofr[h].Ratio(),1);
	for(auto h=0;h<gofr.HisX;h++) {
		auto dd=dx*h;
		cout << h*dx<< "  " <<dd*dd*gofr[h].Ratio() <<endl;
	}

}

void Saxs::readContrast(std::ifstream & in){
	vector<double> x,y;
	double CNtot;
	for(string str;getline(in,str);){
		stringstream sst{str};
		if(str[0] == '#') continue;
		double x0,y0;
		sst>>x0>>y0;
		x.push_back(x0);
		y.push_back(y0);
	}
	double dq=x[1]-x[0];
	double qcut=x[x.size()-1]+dq;
	SaxsHistogram rdf(dq,qcut,unitsQ);
	rdf.setLabel("Contrast");
	for(auto h=0;h<rdf.HisX;h++)
		rdf[h]+=hist1D(y[h],1);
	SContrast=rdf;
};

void Saxs::SetupQdf(){
	if(noSplineOut)
		qdfx=new SaxsHistogram;
	else
		qdfx=new SaxsHistogramSpline;
}

void Saxs::Setup(const vector<string> & at, bool mySans,bool myDens){
	bSans=mySans;
	bElDens=myDens;
	vector<string> at_loc=at;
	sort(at_loc.begin(),at_loc.end());
	auto it=unique(at_loc.begin(),at_loc.end());
	at_loc.resize(distance(at_loc.begin(),it));
	ScatteringFactors * tmp;
	cout << " ---------- Atoms Count ---------- \n\n";
	if(bSans)
		tmp=new ScatteringFactorsN;
	else
		tmp=new ScatteringFactors;
	for(auto atv: at_loc){
		Sfacts[atv]=(*tmp)(atv);
		iSfacts[atv]=opgather<string>{atv}(at);
		cout << " -- Atoms  " <<atv << "    Count =  " << iSfacts[atv].size()<<endl;
	}
	cout << "\n"<<endl;
}
void Saxs::Setup(const vector<int> & Lst,const vector<string> & at1,bool mySans,bool myDens){
	bSans=mySans;
	bElDens=myDens;
	Sfacts.clear();
	iSfacts.clear();
	MyFs.clear();
	MyNas.clear();
	vector<string> at(at1.size(),"0");
	for(auto o1=0;o1<Lst.size();o1++){
		auto o=Lst[o1];
		at[o]=at1[o];
	}
	vector<string> at_loc=at;

	sort(at_loc.begin(),at_loc.end());
	auto it=unique(at_loc.begin(),at_loc.end());
	at_loc.resize(distance(at_loc.begin(),it));
	ScatteringFactors * tmp;
	if(bSans)
		tmp=new ScatteringFactorsN;
	else
		tmp=new ScatteringFactors;

	cout << " ---------- Atoms Count ---------- \n\n";
	for(auto atv: at_loc){
		if(atv == "0") continue;
		Sfacts[atv]=(*tmp)(atv);
		iSfacts[atv]=opgather<string>{atv}(at);
		cout << " -- Atoms\t" <<atv << "    Count =\t" << iSfacts[atv].size()<<endl;
	}
	Ntot=0;
	for(auto it=Sfacts.begin();it!=Sfacts.end();it++){
		auto Na=iSfacts[it->first].size();
		MyFs.push_back((*tmp)(it->first));
		MyNas.push_back(Na);
		Ntot+=Na;
	}
	cout << "\n"<<endl;
	Lst_i=Lst;
	at1_i=at1;
}
void Saxs::Modulus(array3<Complex> & ro_k,array3<Complex> & ro_k1){

#pragma omp parallel for collapse(3)
	for(auto i=0;i<nx;i++){
		for(auto j=0;j<ny;j++){
			for(auto k=0;k<nzp;k++){
				ro_k[i][j][k]*=conj(ro_k1[i][j][k]);
			}
		}
	}

}


SaxsHistogram & SaxsHistogram::operator^=(const SaxsHistogram & z){
	dx=z.dx;
	cut=z.cut;
	Label=z.Label;
	HisX=z.HisX;
	try{
		if(hist.size() == 0){
			hist=vector<hist1D>(z.hist.size(),hist1D());
			} else if(hist.size() != z.hist.size()) throw string("The hist1D dimension are not identical ");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}

	for(size_t o=0; o < hist.size(); o++) {
		double tmp=z.hist[o].nhis?z.hist[o].hisX/static_cast<double>(z.hist[o].nhis):0.0;
		hist[o].hisX+=tmp;
		hist[o].nhis=1;
	}

	return *this;
}
void Saxs::WriteI_r(std::ostream & fout){
// Everything is in Angstroems

	double dx=10.0*MCO[XX][XX]/static_cast<double>(nx);
	double dy=10.0*MCO[YY][YY]/static_cast<double>(ny);
	double dz=10.0*MCO[ZZ][ZZ]/static_cast<double>(nz);
	float dvol=static_cast<float>(dx*dy*dz);
	I_r/=dvol;
// Average I_r over M bons on both sides
	int nx0=static_cast<int>(nx);
	int ny0=static_cast<int>(ny);
	int nz0=static_cast<int>(nz);
	auto Idx=[&nx0](const int mx){return mx<0?mx+nx0:(mx<nx0?mx:mx-nx0);};
	auto Idy=[&ny0](const int my){return my<0?my+ny0:(my<ny0?my:my-ny0);};
	auto Idz=[&nz0](const int mz){return mz<0?mz+nz0:(mz<nz0?mz:mz-nz0);};

	int M=DensAvg;
	double Norm{static_cast<double>((2*M+1)*(2*M+1)*(2*M+1))};
	array3<double> I_r0(nx,ny,nz);
	array3<Complex> ro_k(nx,ny,nzp,align);
	rcfft3d Forward3(nz,I_r,ro_k);
	for(size_t o{0};o<nx;o++)
		for(size_t p{0};p<ny;p++)
			for(size_t q{0};q<nz;q++){
				double Val{0};
				int mj{0};
				for(int n{-M};n<=M;n++){
					auto ox = Idx(o + n);
					for(int m{-M};m<=M;m++){
						auto oy = Idy(p + m);
						for(int l{-M};l<=M;l++){
							auto oz = Idz(q + l);
							Val+=I_r[ox][oy][oz];
							mj++;
						}
					}
				}
				I_r0[o][p][q]=Val/Norm;
			}
// Get the average density at border

	M=4;


	double ISum0{0},ISum1{0};
	for(size_t o{0};o<nx;o++)
		for(size_t p{0};p<ny;p++)
			for(size_t q{0};q<nz;q++){
				ISum0+=I_r0[o][p][q];
			}
	for(int o{M};o<nx-M;o++)
		for(int p{M};p<ny-M;p++)
			for(int q{M};q<nz-M;q++){
				ISum1+=I_r0[o][p][q];
			}
	double Bins=static_cast<double>((nx*ny*nz)-(nx-2*M)*(ny-2*M)*(nz-2*M));
	double AvgDens=(ISum0-ISum1)/Bins;

	int mm0{0},mm1{0};
	for(size_t o{0};o<nx;o++){
		for(size_t p{0};p<ny;p++){
			for(size_t q{0};q<nz;q++){
				I_r0[o][p][q]-=AvgDens;
			}
		}
	}
	int N_t{0};
	for(auto it=Sfacts.begin();it!=Sfacts.end();it++){
		auto Na=iSfacts[it->first].size();
		N_t+=Na;
	}
	I_r=I_r0;
	I_r*=dvol;
	ro_k=(0.0,0.0);
	Forward3.fft(I_r,ro_k);
	Modulus(ro_k,ro_k);
	Complex iD{1.0/N_t,0.0};
	ro_k*=iD;
	I_k=ro_k;
	count=1;
	cout << "\n ------> Average electron density at border is "<< AvgDens << " <-------\n"<<endl;

#ifdef HAVE_VTK
	float x[3], a, v[3], rMin = 0.5, rMax = 1., deltaRad, deltaZ;
	int dims[3]={static_cast<int>(nx),static_cast<int>(ny),static_cast<int>(nz)};
	  // Create the structured grid.
	vtkStructuredGrid *sgrid = vtkStructuredGrid::New();
	sgrid->SetDimensions(dims);

	  // Create points.
	vtkPoints *points = vtkPoints::New();
	points->Allocate(dims[0]*dims[1]*dims[2]);

	  // Create scalars. No need to allocate size, as the array will expand automatically as we add new
	  // values with InsertNextValue.
	vtkFloatArray *scalars = vtkFloatArray::New();
	scalars->Allocate(dims[0]*dims[1]*dims[2]);
	scalars->SetName("density");
	for (size_t k{0}; k < nz; k++) {
		x[ZZ]=dz*k;
		int kOffset = k*nx*ny;
	    for (size_t j{0}; j<ny; j++) {
	    	int jOffset = j * nx;
			x[YY]=dy*j;
	    	for (size_t i{0}; i<nx; i++) {
	    		int offset = i + jOffset + kOffset;
	    		x[XX]=dx*i;

	    		float Val=I_r0[i][j][k];
	    		points->InsertPoint(offset,x);
	    		scalars->InsertValue(offset,Val);
	    	}
	    }
	}
	sgrid->SetPoints(points);
	points->Delete();
	sgrid->GetPointData()->SetScalars(scalars);
	scalars->Delete();
	vtkXMLStructuredGridWriter *Writer = vtkXMLStructuredGridWriter::New();
	Writer->SetInputData(sgrid);
	Writer->SetFileName(fileDens.c_str());
	Writer->Write();

	// Delete objects.
	Writer->Delete();
	sgrid->Delete();

#else
	for(size_t i{0};i<nx;i++)
		for(size_t j{0};j<ny;j++)
			for(size_t k{0};k<nz;k++){
				fout << dx*i<< ", " << dy*j<< ", " << dz*k<< ", " << I_r0[i][j][k] << endl;
			}
#endif
try{
	if(fileDens.find_first_of(".") != std::string::npos){
		fileDens=fileDens.substr(0,fileDens.find_first_of(".")+1)+"qout";
	} else{
		throw string("Output file with no extension. This is weired!!");
	}
}catch(const string & s){
	cout << s <<endl;
	exit(1);
}

	ofstream fout_n(fileDens.c_str(),ios::out);
	if(qdfx == nullptr){
		if(noSplineOut)
			qdfx=new SaxsHistogram;
		else{
			qdfx=new SaxsHistogramSpline;
		}
		*qdfx=__Qhistogram();
	}
	fout_n <<*qdfx ;

}
void Saxs::WriteIt(std::ostream & fout){
	try{
		if(!bAvg) throw string("\nCannot write to file without averaging.");
	} catch(const string & s){cout << s<<endl;Finale::Finalize::Final();}
	if(bElDens){
		WriteI_r(fout);
		return;
	}
	if(qdfx == nullptr){
		if(noSplineOut)
			qdfx=new SaxsHistogram;
		else{
			qdfx=new SaxsHistogramSpline;
		}
		*qdfx=__Qhistogram();
	}
	fout <<*qdfx ;
}

Saxs::~Saxs() {
	// TODO Auto-generated destructor stub
	delete qdfx;
}
std::ostream & operator<<(std::ostream & fout, Saxs & y){
	y.WriteIt(fout);
	return fout;
}
Enums::Padding Saxs::pickPadding(const vector<int> & y) const{
	for(auto op: y) if(op != 0) return Enums::avgDensity;
	return Enums::zero;
}
bool Saxs::HaveSolute(vector<int> & y) const{
	for(auto op: y) if(!op) return true;
	return false;
}

template  <>
void Saxs::Compute<Enums::SQ>(RhoSaxs * x,MAtoms * y){
	ComputeSq(x,y);
}
template  <>
void Saxs::Compute<Enums::SAXS>(RhoSaxs * x,MAtoms * y){
	ComputeSAXS(x,y);
}
template  <>
void Saxs::Compute<Enums::SANS>(RhoSaxs * x,MAtoms * y){
	ComputeSANS(x,y);
}
template  <>
void Saxs::Compute<Enums::ELDENS>(RhoSaxs * x,MAtoms * y){
	ComputeDENS(x,y);
}
void Saxs::ComputeSANS(RhoSaxs * Rho_ex,const MAtoms * y){
	AtomsD * x0{new MAtoms(*y)};
	if(SuperCell0 >1 )
		if(!x0->CenterAtoms()) return;

/*
 * Pick padding. Padding is either average padding from box edges or custom for all systems except
 * for solute-only systems. In that case only zero padding is considered
 */
	Enums::Padding myPadding=Rho_ex->WhichPadding();
	void (RhoSaxs::*Padding)(RhoSaxs *,int,int,int){nullptr};
/*
 * Pick a pointer-to-member function to apply padding according to user's choice
 */
	switch(myPadding){
	case zero:
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
	size_t nnr=x0->getNR();
	vector<vector<int> > g=y->getSaxsSolute();
	int tot{0};
	for(auto op: g)for(auto oq:op)tot++;
	Nsolute=tot;
	RhoSaxs & Rho_e=*Rho_ex;
	array3<Complex> ro_k(nx,ny,nzp,align);
	array3<Complex> ro_ktot(nx,ny,nzp,align);
	RhoSaxs Rho_alt(nx,ny,nz);
	array4<double> & ro_r=Rho_alt;
	crfft3d Backward3(nx,ro_k,ro_r[0]);
	rcfft3d Forward3(nx,ro_r[0],ro_k);
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
	ro_ktot={0.0,0.0};
	int N_t=0;
	Nx=Rho_ex->getnnx();
	Ny=Rho_ex->getnny();
	Nz=Rho_ex->getnnz();
	for(auto o=0;o<DIM;o++) {co[o][o]=SuperCell;oc[o][o]=1.0/SuperCell;}
	map<const string,ScatteringFactorsN::opsfact>::iterator it;
	this->__shift(x0);
//	for(it=Sfacts.begin();it!=Sfacts.end();it++){
//		cout << it->first << " " << it->second() <<endl;
//		auto ff=it->second;
//		Rho_e=0.0;
//		Rho_e.Density(order,x0,iSfacts[it->first],it->first);
//		(Rho_alt.*Padding)(Rho_ex,Nx,Ny,Nz);
//		int nn{0};
//		ro_k=Complex{0.0,0.0};
//		Forward3.fft(ro_r[0],ro_k);
//		auto Na=iSfacts[it->first].size();
//		N_t+=Na;
//
//		double fq=ff();
//#pragma omp parallel for
//		for(auto i=0;i<nx;i++){
//			for(auto j=0;j<ny;j++){
//				for(auto k=0;k<nzp;k++){
//					ro_ktot[i][j][k]+=fq*ro_k[i][j][k];
//				}
//			}
//		}
//	}
	Rho_e=0.0;
	for(it=Sfacts.begin();it!=Sfacts.end();it++){
		auto Na=iSfacts[it->first].size();
		vector<double> wei(Na,it->second());
		Rho_e.Density(order,x0,iSfacts[it->first],it->first,wei);
		N_t+=Na;
	}
	(Rho_alt.*Padding)(Rho_ex,Nx,Ny,Nz);
	int nn{0};
	ro_k=Complex{0.0,0.0};
	Forward3.fft(ro_r[0],ro_k);

#pragma omp parallel for
	for(auto i=0;i<nx;i++){
		for(auto j=0;j<ny;j++){
			for(auto k=0;k<nzp;k++){
				ro_ktot[i][j][k]+=ro_k[i][j][k];
			}
		}
	}

	Modulus(ro_ktot,ro_ktot);
	Complex iD{1.0/N_t,0.0};
	ro_ktot*=iD;
	I_k+=ro_ktot;
	count++;
	delete x0;
}
void Saxs::ComputeSq(RhoSaxs * Rho_ex,const MAtoms * y){
/*
 * Pick padding. Padding is either average padding from box edges or custom for all systems except
 * for solute-only systems. In that case only zero padding is considered
 */
	Enums::Padding myPadding=Rho_ex->WhichPadding();
	void (RhoSaxs::*Padding)(RhoSaxs *,int,int,int){nullptr};
	Padding=&RhoSaxs::doPadding<zero>;

	AtomsD x0;
	Mt=y->getMt();
	Nx=Rho_ex->getnnx();
	Ny=Rho_ex->getnny();
	Nz=Rho_ex->getnnz();

	Matrix oc=Mt.getOC();
	Matrix co=Mt.getCO();
	MCO+=co;
	MOC+=oc;
	Vol+=Mt.getVol();
	const vector<vector<int> > & mCluster=y->getCluster();
	const vector<vector<int> > & mAtoms=y->getAtoms();
	double nato{0};
	vector<double> wei(mCluster.size(),0.0);
	for(size_t o=0;o<mCluster.size();o++){
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			nato+=(double) mAtoms[n].size();
			wei[o]+=(double) mAtoms[n].size();
		}
	}
	for(size_t o=0;o<mCluster.size();o++){
		wei[o]*=(double) mCluster.size()/nato;
	}
	double tot{0.0};
	for(size_t o=0;o<mCluster.size();o++){
		tot+=wei[o];
	}
 	CenterMassD & xcm=y->getCM();
	size_t nr=xcm.Size();

	x0.setDim(nr);
	rvec * x=new rvec[nr];
	for(auto o=0;o<nr;o++)
		for(auto p=0;p<DIM;p++){
			x[o][XX]=co[XX][XX]*xcm[o][XX]+co[XX][YY]*xcm[o][YY]+co[XX][ZZ]*xcm[o][ZZ];
			x[o][YY]=co[YY][XX]*xcm[o][XX]+co[YY][YY]*xcm[o][YY]+co[YY][ZZ]*xcm[o][ZZ];
			x[o][ZZ]=co[ZZ][XX]*xcm[o][XX]+co[ZZ][YY]*xcm[o][YY]+co[ZZ][ZZ]*xcm[o][ZZ];
		}
	x0.setCoord(y->getMt(),x);
	RhoSaxs & Rho_e=*Rho_ex;
	array3<Complex> ro_ktot(nx,ny,nzp,align);
	RhoSaxs Rho_alt(nx,ny,nz);
	array4<double> & ro_r=Rho_alt;

	rcfft3d Forward3(nx,ro_r[0],ro_ktot);
	Complex vt,vt2,v0;
	ro_ktot={0.0,0.0};
	int N_t=y->getCM().Size();

	vector<size_t> indx(N_t);
	int n{0};
	std::generate(indx.begin(),indx.end(),[&n](){return n++;});
	Rho_e=0.0;
	Rho_e.Density(order,&x0,indx,"DUM");
	(Rho_alt.*Padding)(Rho_ex,Nx,Ny,Nz);

	//	Rho_alt=Rho_e;

	ro_ktot=Complex{0.0,0.0};
	Forward3.fft(ro_r[0],ro_ktot);

	Modulus(ro_ktot,ro_ktot);
	Complex iD{1.0/N_t,0.0};
	ro_ktot*=iD;
	I_k+=ro_ktot;
	
	count++;
	Ntot=1;
	delete [] x;
}


void Saxs::ComputeSAXS(RhoSaxs * Rho_ex,const MAtoms * y,bool bDens){
	AtomsD * x0{new MAtoms(*y)};
	if(SuperCell0 >1 )
		if(!x0->CenterAtoms()) return;

/*
 * Pick padding. Padding is either average padding from box edges or custom for all systems except
 * for solute-only systems. In that case only zero padding is considered
 */
	Enums::Padding myPadding=Rho_ex->WhichPadding();
	void (RhoSaxs::*Padding)(RhoSaxs *,int,int,int){nullptr};
/*
 * Pick a pointer-to-member function to apply padding according to user's choice
 */
	switch(myPadding){
	case zero:
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
	size_t nnr=x0->getNR();
	vector<vector<int> > g=y->getSaxsSolute();
	int tot{0};
	for(auto op: g)for(auto oq:op)tot++;
	Nsolute=tot;
	RhoSaxs & Rho_e=*Rho_ex;
	array3<Complex> ro_k(nx,ny,nzp,align);
	array3<Complex> ro_ktot(nx,ny,nzp,align);
	RhoSaxs Rho_alt(nx,ny,nz);
	array4<double> & ro_r=Rho_alt;
	crfft3d Backward3(nx,ro_k,ro_r[0]);
	rcfft3d Forward3(nx,ro_r[0],ro_k);
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
	ro_ktot={0.0,0.0};
	int N_t=0;
	Nx=Rho_ex->getnnx();
	Ny=Rho_ex->getnny();
	Nz=Rho_ex->getnnz();
	for(auto o=0;o<DIM;o++) {co[o][o]=SuperCell;oc[o][o]=1.0/SuperCell;}
	map<const string,ScatteringFactors::opsfact>::iterator it;
	this->__shift(x0);
	for(it=Sfacts.begin();it!=Sfacts.end();it++){
		Rho_e=0.0;
		Rho_e.Density(order,x0,iSfacts[it->first],it->first);
		(Rho_alt.*Padding)(Rho_ex,Nx,Ny,Nz);
		ro_k=Complex{0.0,0.0};
		auto ff=it->second;
		auto Na=iSfacts[it->first].size();
		N_t+=Na;

		Forward3.fft(ro_r[0],ro_k);

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
	delete x0;
}
void Saxs::ComputeDENS(RhoSaxs * Rho_ex,MAtoms * y){

	AtomsD * x0{new MAtoms(*y)};


	//SuperCell0=1;
	//SuperCell=1;

	size_t nnr=x0->getNR();
	vector<vector<int> > g=y->getSaxsSolute();
	int tot{0};
	for(auto op: g)for(auto oq:op)tot++;
	Nsolute=tot;
	RhoSaxs & Rho_e=*Rho_ex;
	array3<Complex> ro_k(nx,ny,nzp,align);
	array3<Complex> ro_ktot(nx,ny,nzp,align);
	array4<double> & ro_r=Rho_e;
	crfft3d Backward3(nx,ro_k,ro_r[0]);
	rcfft3d Forward3(nx,ro_r[0],ro_k);

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
	ro_ktot={0.0,0.0};
	int N_t=0;
	Nx=Rho_ex->getnnx();
	Ny=Rho_ex->getnny();
	Nz=Rho_ex->getnnz();

	map<const string,ScatteringFactors::opsfact>::iterator it;
	this->__shift(x0);
	if(DensPick["R"]){

		for(it=Sfacts.begin();it!=Sfacts.end();it++){
			Rho_e=0.0;
			Rho_e.Density(order,x0,iSfacts[it->first],it->first);
			auto ff=it->second;
			auto myFF=ff(0.0);
			for(size_t i{0};i<nx;i++)
				for(size_t j{0};j<ny;j++)
					for(size_t k{0};k<nz;k++){
						I_r[i][j][k]+=myFF*Rho_e[0][i][j][k];
					}
		}
	}else{
		ro_ktot={0.0,0.0};
		for(it=Sfacts.begin();it!=Sfacts.end();it++){
			Rho_e=0.0;
			Rho_e.Density(order,x0,iSfacts[it->first],it->first);
			ro_k=Complex{0.0,0.0};
			auto ff=it->second;
			Forward3.fft(ro_r[0],ro_k);

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
		Rho_e=0.0;
		Backward3.fftNormalized(ro_ktot,ro_r[0]);
		for(size_t i{0};i<nx;i++)
			for(size_t j{0};j<ny;j++)
				for(size_t k{0};k<nz;k++){
					I_r[i][j][k]+=ro_r[0][i][j][k];
				}


	}
	count++;
	delete x0;
}

void Saxs::Clear(){
	count=0;
	I_k={0.0,0.0};
	I_r=0.0;
	MCO=0.0;
	MOC=0.0;
}
void Saxs::bReadIq(ifstream & fin){
	size_t mx,my,mz,mzp;
	Matrix oc{0.0},co{0.0};
	fin.seekg (0, fin.end);
	int length = fin.tellg();
	fin.seekg (0, fin.beg);
	try{
	  if(!length) throw string("\n Something is wrong here!! File length is zero. \n");
	}catch(const string & s){
	  cout << s <<endl;
	  Finale::Finalize::Final();
	}

	fin.read(as_byte(Ntot_in),sizeof(Ntot_in));
	Ntot=Ntot_in;
	fin.read(as_byte(Nsolute),sizeof(Nsolute));
	fin.read(as_byte(SuperCell),sizeof(SuperCell));
	fin.read(as_byte(co),sizeof(co));
	fin.read(as_byte(oc),sizeof(oc));
	fin.read(as_byte(mx),sizeof(mx));
	fin.read(as_byte(my),sizeof(my));
	fin.read(as_byte(mz),sizeof(mz));
	fin.read(as_byte(mzp),sizeof(mzp));
	this->Allocate(mx,my,mz);
	fin.read(as_byte(I_k[0][0][0]),sizeof(I_k[0][0][0])*nx*ny*nzp);
	MCO=co;
	MOC=oc;
	bAvg=true; // no need to average here
}
void Saxs::bPrintIq(ofstream & fout){
	Matrix oc{MOC},co{MCO};
	fout.write(as_byte(Ntot),sizeof(Ntot));
	fout.write(as_byte(Nsolute),sizeof(Nsolute));
	fout.write(as_byte(SuperCell),sizeof(SuperCell));
	fout.write(as_byte(co),sizeof(co));
	fout.write(as_byte(oc),sizeof(oc));
	fout.write(as_byte(nx),sizeof(nx));
	fout.write(as_byte(ny),sizeof(ny));
	fout.write(as_byte(nz),sizeof(nz));
	fout.write(as_byte(nzp),sizeof(nzp));
	fout.write(as_byte(I_k[0][0][0]),sizeof(I_k[0][0][0])*nx*ny*nzp);
}
