/*
 * Atoms.cpp
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */
#include "Atoms.h"
#include "jacobi.h"
template <typename T>
int Atoms<T>::calls=0;
template <typename T>
int Atoms<T>::cPrint_calls=0;
template <typename T>
int Atoms<T>::cPDBavg_calls=0;



template <typename T>
vector<string> * Atoms<T>::ResList0=nullptr;

template <typename T>
vector<int> * Atoms<T>::ResIndx0=nullptr;

template <typename T>
map<string,double> Atoms<T>::MMass={{"H",1.008},{"C",12.011},{"N",14.0067},{"O",15.9994},{"S",32.06},{"P",30.974},
		{"FE",55.847},{"CA",40.080},{"ZN",65.37},{"NA",22.98977},{"MG",24.305}};


template <typename T>
map<string,double> Atoms<T>::MMassNCH={{"H",0.0},{"C",0.0},{"N",14.0067},{"O",15.9994},{"S",32.06},{"P",30.974},
		{"FE",55.847},{"CA",40.080},{"ZN",65.37},{"NA",22.98977},{"MG",24.305}};

template <typename T>
string Atoms<T>::MIons="CA NA SOD CAL ZN";

template <typename T>
string Atoms<T>::Mdetg="AOT SDS LEA LDA DOP";

template <typename T>
class Jacob{
	double ** iner;
	double ** ei;
	double * Im;
	static Dvect Im_s;
	void  eigsort(DDvect<T>  & d, MMatrix<T>  * v=NULL){
		int k;
		int n=DIM;
		for (int i=0;i<n-1;i++) {
			double p=d[k=i];
			for (int j=i;j<n;j++)
				if (d[j] >= p) p=d[k=j];
			if (k != i) {
				d[k]=d[i];
				d[i]=p;
				if (v != NULL){
					for (int j=0;j<n;j++) {
						p=(*v)[j][i];
						(*v)[j][i]=(*v)[j][k];
						(*v)[j][k]=p;
					}
				}
			}
		}
	}
public:
	Jacob():  iner(NULL), ei(NULL),Im(NULL){};
	~Jacob(){
		for(int o=0;o<DIM;o++) delete [] iner[o];
		delete [] iner;
		for(int o=0;o<DIM;o++) delete [] ei[o];
		delete [] ei;
		delete [] Im;
	}
	void operator()(MMatrix<T> & iner0, DDvect<T> & Im0, MMatrix<T> & ei0, int & nrot){
		if(!iner) {
			iner=new double * [DIM];
			for(int o=0;o<DIM;o++) iner[o]=new double [DIM];
		}
		if(!ei) {
			ei=new double * [DIM];
			for(int o=0;o<DIM;o++) ei[o]=new double [DIM];
		}
		if(!Im) Im=new double [DIM];
		for(int o=0;o<DIM;o++){
			for(int p=0;p<DIM;p++){
				iner[o][p]=iner0[o][p];
			}
		}
		jacobi(iner,DIM,Im,ei,&nrot);

		for(int o=0;o<DIM;o++){
			Im0[o]=Im[o];
			for(int p=0;p<DIM;p++){
				ei0[o][p]=ei[o][p];
			}
		}

		eigsort(Im0,&ei0);
	}

};

template <typename T>
class PBCvect{
	using Dvect=DDvect<T>;
	static vector<Dvect> * nxyz;
	PBCvect(){};
	static void __compVect(){
		const int M{3};
		nxyz=new vector<Dvect>;
		for(int o=0;o<M;o++){
			double nx=static_cast<double>(o-1);
			for(int p=0;p<M;p++){
				double ny=static_cast<double>(p-1);
				for(int q=0;q<M;q++){
					double nz=static_cast<double>(q-1);
					nxyz->push_back(Dvect(nx,ny,nz));
				}
			}
		}

	};
public:
	static vector<Dvect> & getVec(){__compVect();return *nxyz;}
};

template <typename T>
vector<DDvect<T>> * PBCvect<T>::nxyz=nullptr;

template <typename T>
class CellComp{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	Dvect Xd;
	Matrix co;
public:
	CellComp(Dvect& x, Dvect & y, Matrix & c): Xd(x-y), co(c){};
	CellComp(Dvect& x, Matrix & c): Xd(x), co(c){};
	bool operator()(const Dvect & x, const Dvect & y){
		Dvect x1=Xd+x;
		Dvect x2=Xd+y;
		Dvect xc1=co*x1;
		Dvect xc2=co*x2;
		return xc1.Norm() < xc2.Norm();
	}
};
template <typename T>
class RgComp{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	Matrix co;
	vector<Dvect> * X{nullptr};
	Dvect Xd{0};
	Dvect Rcm{0};
	T Rg{0};
	T fact{0};
	T myRg{0};
public:
	RgComp(vector<Dvect> & x, Matrix & c): X(&x), co(c){
		fact=1.0/(T) X->size();
		CompRg();
		};
	void setVect(Dvect x){Xd=x;}
	void CompRg(){
		vector<Dvect> & X0=*X;
		Rg=0;
		Rcm=Dvect{0};
		for(size_t o{0};o<X0.size();o++){
			Dvect x0=co*X0[o];
			Rg+=x0*x0;
			Rcm+=x0;
		}
		Rg*=fact;
		Rcm*=fact;
		myRg=Rg-Rcm*Rcm;
	}
	T getRg(){return myRg;}
	bool operator()(const Dvect & x, const Dvect & y){
		T Rg1{Rg};
		T Rg2{Rg};
		T myRg1{0};
		T myRg2{0};
		Dvect Rcm1{Rcm};
		Dvect Rcm2{Rcm};

		Dvect xc=co*Xd;
		Dvect x1=Xd+x;
		Dvect x2=Xd+y;
		Dvect xc1=co*x1;
		Dvect xc2=co*x2;



		T rg1=xc1*xc1*fact;
		T rg2=xc2*xc2*fact;
		Dvect rcm=xc*fact;
		Rg1+=xc1*xc1*fact - xc*xc*fact;
		Rg2+=xc2*xc2*fact - xc*xc*fact;
		Rcm1+=xc1*fact-xc*fact;
		Rcm2+=xc2*fact-xc*fact;
		myRg1=Rg1-Rcm1*Rcm1;
		myRg2=Rg2-Rcm2*Rcm2;
		return myRg1 < myRg2;
	}
};
template <typename T>
class EnComp{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	Matrix co;
	vector<Dvect> * X{nullptr};
	Dvect Xd{0};
	size_t MMax{0};
public:
	EnComp(vector<Dvect> & x, Matrix & c): X{&x}, co{c}{};
	void setVect(size_t p){MMax=p+1;Xd=(*X)[p];}
	T compEn(){
		vector<Dvect> & X0=*X;
		auto o=MMax-1;
		Dvect xa_o=X0[o];
		Dvect xc_o=co*xa_o;

		T ener_x{0};

		for(size_t p{0};p<MMax-1;p++){
			Dvect xa_p=X0[p];
			Dvect xc_p=co*xa_p;
			Dvect xd=xc_o-xc_p;
			T dist=sqrt(xd[XX]*xd[XX]+xd[YY]*xd[YY]+xd[ZZ]*xd[ZZ]);

			ener_x+=dist;
		}
		return ener_x;

	}
	bool operator()(const Dvect & x, const Dvect & y){
		vector<Dvect> & X0=*X;
		auto o=MMax-1;
		Dvect xa_o1=X0[o]+x;
		Dvect xc_o1=co*xa_o1;
		Dvect xa_o2=X0[o]+y;
		Dvect xc_o2=co*xa_o2;

		T ener_x{0},ener_y{0};

		for(size_t p{0};p<MMax-1;p++){
			Dvect xa_p=X0[p];
			Dvect xc_p=co*xa_p;
			Dvect xd1=xc_o1-xc_p;
			Dvect xd2=xc_o2-xc_p;
			T dist1=sqrt(xd1[XX]*xd1[XX]+xd1[YY]*xd1[YY]+xd1[ZZ]*xd1[ZZ]);
			T dist2=sqrt(xd2[XX]*xd2[XX]+xd2[YY]*xd2[YY]+xd2[ZZ]*xd2[ZZ]);

			ener_x+=dist1;
			ener_y+=dist2;
		}
		return ener_x < ener_y;
	}
};
template <typename T>
Dvect Jacob<T>::Im_s=0.0;

template<typename T>
int    Atoms<T>::step_c=0;

template<typename T>
float   Atoms<T>::prec_c=0.0;

template<typename T>
float Atoms<T>::time_c=0.0;

template<typename T>
float Atoms<T>::getTime(){return time_c;}

template<typename T>
void Atoms<T>::CompCM(){
	CenterMass<T> & R_cm=*R_cmx;
	const vector<vector<int> > & mCluster=Perco->getCluster();
	vector<vector<int> > & mAtoms=Perco->getAtoms();
	vector<Dvect> R_CM=vector<Dvect>(mCluster.size());
	CenterMass<T>::setTime(time_c);
	for(size_t o=0;o<mCluster.size();o++){
		Dvect cm=T{0};
		double tmass=0.0;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int i=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) cm[o1]+=xa[i][o1]*mass[i];
				tmass+=mass[i];
			}
		}
		cm/=tmass;
		R_CM[o]=cm;
	}
	R_cm.setup(R_CM.size()); R_cm(R_CM);
	double a=Mt.getCO()[XX][XX];
	double b=Mt.getCO()[YY][YY];
	double c=Mt.getCO()[ZZ][ZZ];
	R_cm.setAxis(a,b,c);

  //  copy(R_CM.begin(),R_CM.end(),back_inserter(R_cm));
}
template<typename T>
void Atoms<T>::setDim(int n){
	nr=n;
	x=vector<Dvect>(n,T{0});
	xa=vector<Dvect>(n,T{0});
}

template <typename T>
void Atoms<T>::setTopol(Topol_NS::Topol & y){
	const double massH{1.00784};
	const double massC{12.00960};
	const double eps{1.0e-4};
	rd=y.getrd();
	mass=y.getMass();
	massNCH=mass;
	std::transform(massNCH.begin(),massNCH.end(),massNCH.begin(),[=](double x){
		return (abs(x-massH) < eps || abs(x-massC) <eps)?0.0:x;
	});
	PDB=y.gPDB();
};
template <typename T>
void Atoms<T>::initLists(Topol_NS::TopolPDB & y,vector<string> & L){
	PDBs=vector<Topol_NS::PDBdata>(y.Size());

	for(size_t o=0;o<PDBs.size();o++) PDBs[o]=y[o];
	ResList0=new vector<string>(L);
	ResIndx0=new vector<int>;
	int ia=0;
	for(size_t o=0;o<PDBs.size();o++){
		string sub1=PDBs[o].atn;
		string sub2=PDBs[o].resn;
		atres.push_back(sub2);
		atname.push_back(sub1);
		size_t mst=sub1.find_first_not_of("1234")?1:0;

		if(MIons.find(sub2) == string::npos) {
			mass.push_back(MMass[sub1.substr(mst,1)]);
			massNCH.push_back(MMassNCH[sub1.substr(mst,1)]);
		}
		else {
			mass.push_back(MMass[sub1.substr(mst,2)]);
			massNCH.push_back(MMassNCH[sub1.substr(mst,2)]);
		}
		if(!ResList0 || find(ResList0->begin(),ResList0->end(),sub2) != ResList0->end()) {
			ResIndx0->push_back(ia);
		}
		ia++;
	}
}

template <typename T>
template <myAtoms MM>
inline void Atoms<T>::InitSelection(vector<string> & y, Topol_NS::Topol & MyTop)
{
	vector<vector<int> > * mySel;

	try{
		switch(MM){
		case Selection:
			mySel=&SelRes;
			break;
		case Reference:
			mySel=&SaxsSolute;
			break;
		case fftPadding:
			mySel=&myPadding;
			break;
		default:
			throw string(" ");
		}
	}catch(const string & s){cout << s <<endl;exit(1);}

	vector<int> Type=MyTop.getAtomTypeNo();
	try{
		for(size_t o=0;o<y.size();o++){
			if(!MyTop.CheckResidue(y[o])) throw string(" Selected residue " + y[o] + " does not exist in the system. Abort");
		}
	} catch(const string & s){
		cout << s << endl;
		exit(1);
	}
	typedef map<string,vector<vector<int> > >  ResMap;
	typedef ResMap::iterator ResIter;

	ResMap & ResDef=MyTop.getDef();
	string str=std::accumulate(y.begin(),y.end(),string(""),[](string o, string t){return o+" "+t;});
	for(auto o=0;o< MyTop.ResSize();o++){
		if(str.find(MyTop.ResInfo(o)) != string::npos){
			const vector<int> Res=MyTop.nAtomRes(o);
			mySel->push_back(MyTop.nAtomRes(o));
			if(MM == Selection){
				for(auto p=0;p<Res.size();p++){
					TypeNo.push_back(Type[Res[p]]);
				}
			}
		}

	}
}

template <typename T>
bool Atoms<T>::CenterAtoms(){
	vector<bool> bB(nr,false);
	vector<vector<int> > & ind=SaxsSolute;
	double m{0.0};
	for(auto op: ind)
		for(auto qp: op){
			bB[qp]=true;
			m+=1.0;
		}
	if(!m) return true;
	double Mx=(double) m*3;

	Matrix co{Mt.getCO()};
	double Nx{0};

	for(auto o=0;o<nr;o++){
		for(int r=0;r<DIM;r++) {
			float xd=xa[o][r];
			xd=rint(xd-0.5);
			if(bB[o])
				if(xd != 0.0) {
					Nx+=1.0;
				}
		}
	}
	try{
		if(Nx/Mx > 0.03) throw string("Skipping trajectory point.");
	}catch(const string & s){cout << s<<endl;return false;}
	return true;
}

template <typename T>
Atoms<T>::Atoms(const int nind): R_cmx{new CenterMassBW3<T>} {
	nr=nind;
	if(nr) {
		x=vector<Dvect>(nr);
		xa=vector<Dvect>(nr);
	}
}
template <typename T>
Atoms<T>::Atoms(const AtomIndex & id): R_cmx{new CenterMassBW3<T>}{
	nr=id.getN();
	if(nr) {
		x=vector<Dvect>(nr);
		xa=vector<Dvect>(nr);
	}
}


template <typename T>
Atoms<T> & Atoms<T>::operator()(const int natoms){
	nr=natoms;
	if(!x.empty()) x.clear();
	if(!xa.empty()) xa.clear();
	x=vector<Dvect>(nr);
	xa=vector<Dvect>(nr);

	return *this;
};
template <typename T>
void Atoms<T>::CalcGyro(vector<double> & massa,vector<Gyration<T>*> & Rg){
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();
	for(size_t o=0;o<mCluster.size();o++){
		*Rg[o]=0.0;
		Dvect cm,cmG;
		double unit_nmm2=1.0/(unit_nm*unit_nm);

		int ntot=0;
		double tmass=0.0;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int i=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) cm[o1]+=x[i][o1]*massa[i];
				for(int o1=0;o1<DIM;o1++) cmG[o1]+=x[i][o1];
				tmass+=massa[i];
				ntot++;
			}
		}
		cm/=tmass;
		cmG/=static_cast<int> (ntot);
		Dvect x0;
		Dvect Gm,Im,axis;
		Matrix Giner;
		Matrix ei,ei2;
		double MyRg=0;
		int oa=0;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int nn=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) x0[o1]=cmG[o1]-x[nn][o1];
				Giner+=x0%x0*unit_nmm2;
			}
		}
		Giner/=static_cast<double> (ntot);
		int nrot;
		Jacob<T>()(Giner,Gm,ei,nrot);
		for(int o1=0;o1<DIM;o1++) MyRg+=Gm[o1];
		Matrix Inertia;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int nn=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) x0[o1]=cm[o1]-x[nn][o1];

				Inertia[XX][XX]+=(sqr(x0[YY])+sqr(x0[ZZ]))*massa[nn];
				Inertia[YY][YY]+=(sqr(x0[ZZ])+sqr(x0[XX]))*massa[nn];
				Inertia[ZZ][ZZ]+=(sqr(x0[XX])+sqr(x0[YY]))*massa[nn];
				Inertia[XX][YY]-=(x0[XX]*x0[YY])*massa[nn];
				Inertia[XX][ZZ]-=(x0[XX]*x0[ZZ])*massa[nn];
				Inertia[YY][ZZ]-=(x0[YY]*x0[ZZ])*massa[nn];
			}
		}
		Inertia[YY][XX]=Inertia[XX][YY];
		Inertia[ZZ][XX]=Inertia[XX][ZZ];
		Inertia[ZZ][YY]=Inertia[YY][ZZ];
		Jacob<T>()(Inertia,Im,ei2,nrot);
		Im*=unit_nmm2;
		double fact=5.0/tmass/2.0;
		axis[XX]=fact*(Im[YY]+Im[ZZ]-Im[XX]);
		axis[YY]=fact*(Im[XX]+Im[ZZ]-Im[YY]);
		axis[ZZ]=fact*(Im[YY]+Im[XX]-Im[ZZ]);
		MyRg=(axis[XX]+axis[YY]+axis[ZZ])/5.0;
		(*Rg[o])(MyRg,Im,Gm,axis,cm);

	}

}
template <typename T>
template<Enums::myWriteOptions OPT>
void Atoms<T>::Gyro(){
	vector<vector<int> >  mCluster=Perco->getCluster();
	vector<vector<int> >  mAtoms=Perco->getAtoms();
	vector<string>  & Types=Perco->getpResn();
	std::hash<string> str_hash;

	vector<Gyration<T> *> Rg=vector<Gyration<T>*>(mCluster.size());
	switch(OPT){
	case Enums::JSON:
		for(auto & ip: Rg)
			ip=new GyrationJSON<T>();
		break;
	default:
		for(auto & ip: Rg)
			ip=new Gyration<T>();
		break;
	}
	Rg_i.clear();
	Gyration<T>::setTime(time_c);

	switch(OPT){
	case Enums::JSON:
		CalcGyro(mass,Rg);
		for(size_t o{0};o<Rg.size();o++){
			map<string,vector<int>> Tag;
			for(size_t p=0;p<mCluster[o].size();p++){
				int n=mCluster[o][p];
				Tag[Types[n]].push_back(n);
			}
			map<string,size_t> tstr;
			for(auto tags: Tag){
				vector<int> v=tags.second;
				std::sort(v.begin(),v.end());
				std::stringstream ss;
				std::copy(v.begin(),v.end(),std::ostream_iterator<int>( ss," "));
				tstr[tags.first]=str_hash(ss.str());
			}
			double a{Rg[o]->gRadg()};
			Dvect I{Rg[o]->gI()};
			Dvect G{Rg[o]->gG()};
			Dvect axis{Rg[o]->gaxis()};
			Dvect xcm{Rg[o]->gXcm()};
			Rg_i.push_back(new GyrationJSON<T>(a,I,G,axis,xcm,tstr));
		}
		CalcGyro(massNCH,Rg);
		for(size_t o{0};o<Rg.size();o++){
			map<string,vector<int>> Tag;
			for(size_t p=0;p<mCluster[o].size();p++){
				int n=mCluster[o][p];
				Tag[Types[n]].push_back(n);
			}
			map<string,size_t> tstr;
			for(auto tags: Tag){
				vector<int> v=tags.second;
				std::sort(v.begin(),v.end());
				std::stringstream ss;
				std::copy(v.begin(),v.end(),std::ostream_iterator<int>( ss," "));
				tstr[tags.first]=str_hash(ss.str());
			}
			double a{Rg[o]->gRadg()};
			Dvect I{Rg[o]->gI()};
			Dvect G{Rg[o]->gG()};
			Dvect axis{Rg[o]->gaxis()};
			Dvect xcm{Rg[o]->gXcm()};
			Rg_i.push_back(new GyrationJSON<T>(a,I,G,axis,xcm,tstr));
		}
		break;
	default:
		CalcGyro(mass,Rg);
		for(size_t o{0};o<Rg.size();o++){
			Rg_i.push_back(new Gyration<T>(*Rg[o]));
		}
		CalcGyro(massNCH,Rg);
		for(size_t o{0};o<Rg.size();o++)
			Rg_i.push_back(new Gyration<T>(*Rg[o]));
		break;

	}

	Rg_count++;
}

template <typename T>
Atoms<T> & Atoms<T>::operator=(const Atoms<T> & y){
	try{
		if(nr != y.nr) throw "Cannot copy Atoms of different size ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	for(int i=0;i<nr;i++) {
		x[i][0]=y.x[i][0];
		x[i][1]=y.x[i][1];
		x[i][2]=y.x[i][2];
	}
	Mt(y.Mt);
	return *this;
}
template <typename T>
Atoms<T> & Atoms<T>::operator=(const T a){
	for(int i=0;i<nr;i++) {
		x[i][0]=a;
		x[i][1]=a;
		x[i][2]=a;
	}
	return *this;
}
template <typename T>
void Atoms<T>::Rot(const Matrix coR){
	rvec t;
	for(int i=0;i<nr;i++){
		for(int o=0;o<DIM;o++) t[o]=coR[o][0]*x[i][0]+coR[o][1]*x[i][1]+coR[o][2]*x[i][2];
		for(int o=0;o<DIM;o++) x[i][o]=t[o];
	}
}
template <typename T>
void Atoms<T>::setCoord(const Metric<T> & Mt_in, const rvec * x0){
	Mt(Mt_in);
	for(int i=0;i<nr;i++)
		for(int j=0;j<DIM;j++)
			x[i][j]=x0[i][j];
}
template <typename T>
void Atoms<T>::setMT(const Metric<T> & Mt_in){
		Mt(Mt_in);
	}

template <typename T>
void Atoms<T>::doCOtoOC(){
	try{
	if(x.empty()) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	if(xa.empty()) xa=vector<Dvect>(nr);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) xa[i][o]=Mt.getOC()[o][0]*x[i][0]+Mt.getOC()[o][1]*x[i][1]+
			Mt.getOC()[o][2]*x[i][2];
}
template <typename T>
void Atoms<T>::doOCtoCO(){
	try{
	if(xa.empty() || x.empty()) throw "Atoms class atom coordinates undefined";}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) x[i][o]=Mt.getCO()[o][0]*xa[i][0]+Mt.getCO()[o][1]*xa[i][1]+
			Mt.getCO()[o][2]*xa[i][2];
}
template <typename T>
Atoms<T> Atoms<T>::COtoOC(){
	try{
	if(x.empty()) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	Atoms<T> other(nr);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) other.x[i][o]=Mt.getOC()[o][0]*x[i][0]+Mt.getOC()[o][1]*x[i][1]+Mt.getOC()[o][2]*x[i][2];
	return other;
}
template <typename T>
Atoms<T> Atoms<T>::OCtoCO(){
	try{
	if(x.empty()) throw "Atoms class atom coordinates undefined";
	}
	catch(const char *s ){
		std::cout << s << std::endl;
		exit(1);
	}
	Atoms<T> other(nr);
	for(int i=0;i<nr;i++)
		for(int o=0;o<DIM;o++) other.x[i][o]=Mt.getCO()[o][0]*x[i][0]+Mt.getCO()[o][1]*x[i][1]+Mt.getCO()[o][2]*x[i][2];
	return other;
}

template <typename T>
void Atoms<T>::ndxPrint(ostream & fout){
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();
	stringstream ss;
	ss<<time_c;
	fout << string("# Time = "+ ss.str()) <<endl;
	for(size_t r=0;r<mCluster.size();r++){
		fout << "[ Cluster " << r+1 <<" ]"<<endl;
		vector<int> myndx;
		for(size_t s=0;s<mCluster[r].size();s++){
			int t=mCluster[r][s];
			for(size_t w=0;w<mAtoms[t].size();w++){
				string tmp;
				int n=mAtoms[t][w];
				myndx.push_back(n);
			}
		}
		int N{0};
		for( auto idx : myndx){
		  if((N+1)%16)
				fout << idx+1 << " " ;
			else
				if(N) fout <<"\n";

			N++;
		}

		fout << endl;
	}
}
template <typename T>
void Atoms<T>::PrintAll(ostream & fout){
	vector<Dvect> xc=vector<Dvect>(this->nr);

	for(int n=0;n<nr;n++)
		xc[n]=x[n];


	stringstream ss;
	ss<<time_c;
	fout << string("REMARK    Generated by trjProp ") <<endl;
	fout << string("REMARK    SIMULATION TIME = "+ ss.str()) <<endl;

	vector<T> Par=Mt.getParas();
	ss.str(string());
	ss<< setw(9) << setprecision(3) << fixed << Par[0]*10.0;
	ss<< setw(9) << setprecision(3) << fixed << Par[1]*10.0 ;
	ss<< setw(9) << setprecision(3) << fixed << Par[2]*10.0 ;
	ss<< setw(7) << setprecision(2) << fixed << Par[3] ;
	ss<< setw(7) << setprecision(2) << fixed << Par[4] ;
	ss<< setw(7) << setprecision(2) << fixed << Par[5] ;
	fout<< string("CRYST1"+ ss.str() + " P 1           1") <<endl;
	fout << "MODEL "<< fixed << setw(5) << right << cPrint_calls<<endl;
	for(size_t n{0};n<nr;n++){
		string tmp;
		ss.str(string());
		ss <<std::setw(5) << std::left << PDB[n].resn;
		string tmp2=PDB[n].first;
		tmp2.erase(26,1);
		tmp=tmp2.replace(21,5,ss.str());
		ss.str(string());
		ss << std::setw(5) << std::right << PDB[n].res+1;
		tmp.replace(21,5,ss.str());
		tmp.append(" 1.00");
		tmp.append(PDB[n].last);
		ss.str(string());
		ss << fixed << setw(8) << setprecision(3) << right << xc[n][XX]/unit_nm;;
		ss << fixed << setw(8) << setprecision(3) << right << xc[n][YY]/unit_nm;;
		ss << fixed << setw(8) << setprecision(3) << right << xc[n][ZZ]/unit_nm;;
		tmp.replace(30,24,ss.str());
		fout << tmp<<endl;

	}
}
template <typename T>
void Atoms<T>::PDBavg(){
	const vector<vector<int> > & mCluster=Perco->getCluster();
	const vector<vector<int> > & mAtoms=Perco->getAtoms();

	if(cPDBavg_calls == 0){
		size_t Ntot{0};
		for(size_t r=0;r<mCluster.size();r++){
			for(size_t s=0;s<mCluster[r].size();s++){
				int t=mCluster[r][s];
				Ntot+=mAtoms[t].size();
			}
		}
		x_avg=vector<Dvect>(Ntot);
	}

	cPDBavg_calls++;


	vector<Dvect> xc=vector<Dvect>(this->nr);

	size_t N{0};
	for(size_t r=0;r<mCluster.size();r++){
		for(size_t s=0;s<mCluster[r].size();s++){
			int t=mCluster[r][s];
			for(size_t w=0;w<mAtoms[t].size();w++){
				int n=mAtoms[t][w];
				x_avg[N]+=x[n];
				N++;
			}
		}
	}
	nCO+=Mt.getCO();
}
template <typename T>
void Atoms<T>::printPDBavg(ostream & fout){
	const vector<vector<int> > & mCluster=Perco->getCluster();
	const vector<vector<int> > & mAtoms=Perco->getAtoms();
	stringstream ss;
	ss<<time_c;
	fout << string("REMARK    Generated by trjProp ") <<endl;
	fout << string("REMARK    SIMULATION TIME = "+ ss.str()) <<endl;

	vector<T> Par=Mt.getParas();
	ss.str(string());
	ss<< setw(9) << setprecision(3) << fixed << Par[0]*10.0;
	ss<< setw(9) << setprecision(3) << fixed << Par[1]*10.0 ;
	ss<< setw(9) << setprecision(3) << fixed << Par[2]*10.0 ;
	ss<< setw(7) << setprecision(2) << fixed << Par[3] ;
	ss<< setw(7) << setprecision(2) << fixed << Par[4] ;
	ss<< setw(7) << setprecision(2) << fixed << Par[5] ;
	fout<< string("CRYST1"+ ss.str() + " P 1           1") <<endl;

	fout << "MODEL "<< fixed << setw(5) << right << cPrint_calls<<endl;
	size_t N{0};
	for(size_t r=0;r<mCluster.size();r++){
		for(size_t s=0;s<mCluster[r].size();s++){
			int t=mCluster[r][s];
			for(size_t w=0;w<mAtoms[t].size();w++){

				string tmp;
				int n=mAtoms[t][w];
				ss.str(string());
				ss <<std::setw(5) << std::left << PDB[n].resn;
				string tmp2=PDB[n].first;
				tmp2.erase(26,1);
				tmp=tmp2.replace(21,5,ss.str());
				ss.str(string());
				ss << std::setw(5) << std::right << PDB[n].res+1;
				tmp.replace(21,5,ss.str());
				tmp.append(" 1.00");
				tmp.append(PDB[n].last);
				ss.str(string());
				ss << fixed << setw(8) << setprecision(3) << right << x_avg[N][XX]/static_cast<double>(cPDBavg_calls)/unit_nm;
				ss << fixed << setw(8) << setprecision(3) << right << x_avg[N][YY]/static_cast<double>(cPDBavg_calls)/unit_nm;
				ss << fixed << setw(8) << setprecision(3) << right << x_avg[N][ZZ]/static_cast<double>(cPDBavg_calls)/unit_nm;
				tmp.replace(30,24,ss.str());
				fout << tmp<<endl;
				N++;
			}
		}
	}
}

template <typename T>
void Atoms<T>::cPrint(ostream & fout){

	if(bndx) return ndxPrint(fout);

	cPrint_calls++;

	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();


	vector<Dvect> xc=vector<Dvect>(this->nr);

	for(int n=0;n<nr;n++)
		xc[n]=x[n];


	stringstream ss;
	ss<<time_c;
	fout << string("REMARK    Generated by trjProp ") <<endl;
	fout << string("REMARK    SIMULATION TIME = "+ ss.str()) <<endl;

	vector<T> Par=Mt.getParas();
	ss.str(string());
	ss<< setw(9) << setprecision(3) << fixed << Par[0]*10.0;
	ss<< setw(9) << setprecision(3) << fixed << Par[1]*10.0 ;
	ss<< setw(9) << setprecision(3) << fixed << Par[2]*10.0 ;
	ss<< setw(7) << setprecision(2) << fixed << Par[3] ;
	ss<< setw(7) << setprecision(2) << fixed << Par[4] ;
	ss<< setw(7) << setprecision(2) << fixed << Par[5] ;
	fout<< string("CRYST1"+ ss.str() + " P 1           1") <<endl;

	map<int,vector<string> >MapRes;
	fout << "MODEL "<< fixed << setw(5) << right << cPrint_calls<<endl;
	for(size_t r=0;r<mCluster.size();r++){
		for(size_t s=0;s<mCluster[r].size();s++){
			int t=mCluster[r][s];
			for(size_t w=0;w<mAtoms[t].size();w++){

				string tmp;
				int n=mAtoms[t][w];
				ss.str(string());
				ss <<std::setw(5) << std::left << PDB[n].resn;
				string tmp2=PDB[n].first;
				tmp2.erase(26,1);
				tmp=tmp2.replace(21,5,ss.str());
				ss.str(string());
				ss << std::setw(5) << std::right << PDB[n].res+1;
				tmp.replace(21,5,ss.str());
				tmp.append(" 1.00");
				tmp.append(PDB[n].last);
				ss.str(string());
				ss << fixed << setw(8) << setprecision(3) << right << xc[n][XX]/unit_nm;;
				ss << fixed << setw(8) << setprecision(3) << right << xc[n][YY]/unit_nm;;
				ss << fixed << setw(8) << setprecision(3) << right << xc[n][ZZ]/unit_nm;;
				tmp.replace(30,24,ss.str());
				fout << tmp<<endl;
			}
		}

		fout << string("TER") <<endl;
	}
	fout << string("ENDMDL") <<endl;
}

template <typename T>
void Atoms<T>::WriteaStep(FstreamF * fout){
	try{
		throw string("Don't know how to write dcd files. Maybe I should take a class... Abort");
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
}
template <typename T>
void Atoms<T>::WriteaStep(FstreamC * fout){
	XDRFILE * xd=fout->getfin();


	float prec=prec_c;
	matrix box;
	rvec * x0=new rvec[nr];
	for(int o=0;o<nr;o++)
		for(int q=0;q<DIM;q++)
			x0[o][q]=static_cast<float>(x[o][q]);


	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++)
			box[i][j]=Mt.getCO()[i][j];
	try{
		if(write_xtc(xd,nr,step_c,time_c,box,x0,prec)) throw string("Cannot write next frame. Abort");
	}catch(const string & s){
		cout << s << endl;
		exit(1);
	}
	delete [] x0;

}

template <typename T>
void Atoms<T>::ReadaStep(std::ifstream & fin){
	std::cout << " Why am I here?? "<<  std::endl;
	exit(1);
}

template <typename T>
void Atoms<T>::ReadaStep(FstreamC * fin){
	XDRFILE * xd=fin->getfin();
	rvec   *x0;
	matrix box;
	int    step;
	float   prec,time;
	int my_nframe=fin->gFrame();

	FILE * fp=xdrfile_get_fp(xd);
	XDR * myxdr=xdrfile_get_xdr(xd);
	try{
		if(firsttime){
			firsttime=false;
			if(my_nframe != fin->goffStep()){
				if(!(my_nframe<= fin->gFrameNumber()))
					throw string(" Beginning frame is out of range ");
				}
			fin->Rewind();
		}
		x0=new rvec[nr];
		// Apparenly xdr_xtc_seek_frame does not work for the beginning frame!!
		if(my_nframe!=fin->goffStep()) xdr_xtc_seek_frame(my_nframe,fp,myxdr,nr);
		if(read_xtc(xd,nr,&step,&time,box,x0,&prec))
			throw string("Cannot read next frame!");
	}
	catch(const string & s){
		cout << s << endl;
		exit(-1);
	}
	prec_c=prec;
	time_c=time;
	step_c=step;
	MMatrix<T> box1;
	for(auto o=0;o<DIM;o++)
		for(auto p=0;p<DIM;p++)
			box1[o][p]=static_cast<T> (box[o][p]);
	Metric<T> Met(box1);
	setCoord(Met,x0);
	doCOtoOC();
	fin->nextFrame();
	delete []x0;
}

template <typename T>
void Atoms<T>::ReadaStep(FstreamF * fin){
	ReadaStepF(fin->getfin());
}

template <typename T>
void Atoms<T>::pdb(const vector<string> & data){

	vector<Dvect> xc;
	double aa=-1.0,bb=-1.0,cc=-1.0,alpha=90.0,beta=90.0,gamma=90.0;
	for(size_t i=0;i<data.size();i++){
		if(data[i].find("CRYST1") == 0){

			istringstream ss(data[i].substr(6,49));
			ss>>aa>>bb>>cc>>alpha>>beta>>gamma;
			aa*=unit_nm;
			bb*=unit_nm;
			cc*=unit_nm;
		}
	}
	try{
	if(aa < 0 || bb < 0 || cc <0) throw string("\nWarning: PDB file does not have CRYST1 keyword. ")+
			string("If you are doing calculations\n     on this structure, box is not set.\n");
	} catch(const string & s){cout << s<<endl;};
	Metric<T> Met(brot<T>()(aa,bb,cc,alpha,beta,gamma));
	for(size_t i=0;i<data.size();i++){
		if(data[i].find("ATOM") == 0 || data[i].find("HETATM") == 0 ) {
			Dvect tmp;
			std::stringstream ss(data[i].substr(30,24));
			ss>> tmp;
			xc.push_back(tmp);
		}
	}
	try{
		if(this->nr != int(xc.size())) throw string("Something wrong: No. atoms do not match");
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}

	rvec * x0=new rvec[this->nr];
	for(int n=0;n<this->nr;n++){
		for(int o=0;o<DIM;o++) x0[n][o]=xc[n][o]*unit_nm;
	}
	this->setCoord(Met,x0);
	this->doCOtoOC();
	delete [] x0;
}


template <typename T>
void Atoms<T>::ReadaStepF(std::ifstream & fin){
	double table[6];
	fin.seekg(FORTRANBYTES,ios::cur);
	fin.read(reinterpret_cast<char *> (&table),sizeof(table));
	fin.seekg(FORTRANBYTES,ios::cur);
	table[0]*=unit_nm;
	table[2]*=unit_nm;
	table[5]*=unit_nm;
	double alpha = 90.0 - asin(table[1]) * 90.0 / M_PI_2;
	double beta = 90.0 - asin(table[3]) * 90.0 / M_PI_2;
	double gamma = 90.0 - asin(table[4]) * 90.0 / M_PI_2;

	Metric<T> Met(brot<T>()(table[0],table[2],table[5],alpha,beta,gamma));
	vector<float> X[DIM];

	for(int i=0;i<DIM;i++){
		X[i].resize(nr,0.0);
		fin.seekg(FORTRANBYTES,ios::cur);
		fin.read(reinterpret_cast<char *> (&X[i][0]),sizeof(X[0][0])*nr);
		fin.seekg(FORTRANBYTES,ios::cur);
	}
	rvec * x=new rvec[nr];
	for(int n=0;n<nr;n++){
		for(int o=0;o<DIM;o++) x[n][o]=X[o][n]*unit_nm;
	}
	setCoord(Met,x);
	doCOtoOC();
	time_c+=1.0;

}
template <typename T>
void Atoms<T>::moveOffset(FstreamC * fin){
	fin->nextFrame();
}
template <typename T>
void Atoms<T>::moveOffset(FstreamF * fin){
	moveOffset(fin->getfin());
}
template <typename T>
void Atoms<T>::moveOffset(std::ifstream & fin){
	const int DOUBLE=8;
	const int FLOAT=4;
	ios::streamoff OFFSET=FORTRANBYTES*2+6*DOUBLE+(FORTRANBYTES*2+nr*FLOAT)*DIM;
	fin.seekg(OFFSET,ios::cur);
}

template <typename T>
template<Enums::myWriteOptions OPT>
void Atoms<T>::SetupPercolate(Topol_NS::Topol & myTop){
	auto & Reference=myTop.gReferenceResidues();

	auto & atmss=myTop.getAtomName();
	auto & Index=myTop.gCIndex();

	vector<string> resn(nr);
	for(int o{0};o<nr;o++)
		resn[o]=myTop.AtomResidue(o);
	vector<vector<int>> MySel;
	for(size_t o{0};o<Reference.size();o++)
		MySel.push_back(Index[Reference[o]]);
	switch(OPT){
	case Enums::JSON:
		this->Perco=new PercolationJSON<T>(MySel,rd,resn,atmss);
		break;
	default:
		this->Perco=new Percolation<T>(MySel,rd,resn,atmss);
		break;
	}
}

template <typename T>
void Atoms<T>::SetupPercolate(){
	string Exclusion{"ISU ISO"};
	vector<string> resn=vector<string>(PDBs.size());
	vector<vector<int> > MySel;
	for(size_t o=0;o<resn.size();o++){resn[o]=PDBs[o].resn;}
	vector<string> atmss=vector<string>(nr);
	for(size_t o=0;o<PDBs.size();o++)
		atmss[o]=PDBs[o].atn;
	for(auto o=0; o<SaxsSolute.size();o++){
		auto i=SaxsSolute[o][0];
		if(Exclusion.find(resn[i]) != string::npos)
			continue;
		MySel.push_back(SaxsSolute[o]);
	}
	Perco=new Percolation<T>(MySel,rd,resn,atmss);
}

template <typename T>
int Atoms<T>::Percolate() {
	try{
		if(!this->Perco) throw string("Should initialize percolation. Abort.");
	} catch(const string & s){
		cout << s <<endl;
		exit(1);
	}
	vector<Dvect> v{x};
	Matrix co(Mt.getCO());
	Matrix oc(Mt.getOC());
	this->Perco->doContacts(v,co,oc);
	int result=this->Perco->gCluster();
	this->Perco->Accumulate();
	return result;
}
template <typename T>
const vector<vector<int> > & Atoms<T>::getCluster() const{
	try{
		if(!this->Perco) throw string("This cannot be called if Percolation is not done! ");
	}catch(const string & s){cout <<  s <<endl;exit(0);}
	return Perco->getCluster();
}

template <typename T>
const vector<vector<int> > & Atoms<T>::getAtoms() const{
	try{
		if(!this->Perco) throw string("This cannot be called if Percolation is not done! ");
	}catch(const string & s){cout <<  s <<endl;exit(0);}
	return Perco->getAtoms();
}

template <typename T>
vector<DDvect<T>> Atoms<T>::getGC(){
	try{
		if(!this->hasPerco()) throw string("Cannot compute clusters center of mass without percolation.");
	}catch(const string & s){
		cout << s<<endl;exit(1);
	}
	const vector<vector<int> > & mCluster=Perco->getCluster();
	vector<vector<int> > & mAtoms=Perco->getAtoms();
	vector<Dvect> R_CM=vector<Dvect>(mCluster.size());
	for(size_t o=0;o<mCluster.size();o++){
		Dvect cm{0};
		T tmass=0.0;
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t q=0;q<mAtoms[n].size();q++){
				int i=mAtoms[n][q];
				for(int o1=0;o1<DIM;o1++) cm[o1]+=xa[i][o1];
				tmass+=1.0;
			}
		}
		cm/=tmass;
		R_CM[o]=cm;
	}
	return R_CM;
}
template<typename T>
void Atoms<T>::__ReconstructOneCluster(vector<bool> & atSolv){
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();
	Matrix co=Mt.getCO();
	Matrix oc=Mt.getOC();
	Dvect xcmC{0};
	size_t u{0};
	for(size_t p=0;p<mCluster[0].size();p++){
		int n=mCluster[0][p];
		for(size_t i=0;i<mAtoms[n].size();i++){
			int ia=mAtoms[n][i];
			xcmC[XX]+=xa[ia][XX];
			xcmC[YY]+=xa[ia][YY];
			xcmC[ZZ]+=xa[ia][ZZ];
			u++;
		}
	}
	xcmC/=static_cast<double>(u);

	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			 xa[ia][o]-=xcmC[o]-HALF;
			 if(atSolv[ia]) {
				 xa[ia][o]-=rint(xa[ia][o]-HALF);
			 }
		}
	}
	//> Obtain new Cartesian coordinates from reduced xa's
	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			x[ia][o]=Mt.getCO()[o][XX]*xa[ia][XX]+Mt.getCO()[o][YY]*xa[ia][YY]+Mt.getCO()[o][ZZ]*xa[ia][ZZ];
		}
	}
	//> Obtain new Cartesian coordinates from reduced xa's
	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			x[ia][o]=Mt.getCO()[o][XX]*xa[ia][XX]+Mt.getCO()[o][YY]*xa[ia][YY]+Mt.getCO()[o][ZZ]*xa[ia][ZZ];
		}
	}
	//this->PrintAll(cout);exit(1);
}
template <typename T>
void Atoms<T>::Reconstruct(Contacts<T> * con0){
	vector<Dvect> nxyz=PBCvect<T>::getVec();
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();

	// Lambda function find minimum Rg for a list of vectors
		// Lambda function find minimum Rg for a list of vectors
	auto cmSweep=[nxyz](vector<Dvect> & x,size_t size,EnComp<T> & enComp){
		bool notOk{true};
		for(size_t p=1;p<size;p++){
			enComp.setVect(p);
			Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),enComp);
			for(int q=0; q < DIM;q++) x[p][q]=x[p][q]+tmp[q];
		}
	};


	Matrix co=Mt.getCO();
	Matrix oc=Mt.getOC();

// All coordinates are in reduced units

	Dvect xcmCell{0}; // Geometric center of the aggregate composed of clusters
	vector<Dvect> xcmC(mCluster.size()); // Geometric center of the clusters
	vector<Dvect> xcm(mAtoms.size(),Dvect{T{0.0}}); // Geometric center of the solute residues
	vector<Dvect> xb(nr,Dvect{T{0.0}}); // Atomic coordinates relative to the residue geometric center

	vector<bool> atSolv(nr,true);

	// find atoms which are solvent

	for(auto o=0;o<mAtoms.size();o++){
		for(auto p=0;p<mAtoms[o].size();p++){
			atSolv[mAtoms[o][p]]=false;
		}
	}

// xb is the new coordinate referred to residue geometric center

	for(size_t o=0;o<mAtoms.size();o++){
		vector<Dvect> xa0(mAtoms[o].size(),Dvect{T{0}});
		xcm[o]=Dvect{0};
		for(size_t p=0;p<mAtoms[o].size();p++){
			int n=mAtoms[o][p];
			xa0[p]=xa[n];
		}
		auto Encmp=EnComp<T>(xa0,co);
		cmSweep(xa0,mAtoms[o].size(),Encmp);

		for(size_t p=0;p<mAtoms[o].size();p++){
			Dvect xx=xa0[p];
			xcm[o]+=xx;
		}

		xcm[o]/=static_cast<double>(mAtoms[o].size());

		for(size_t p=0;p<mAtoms[o].size();p++){
			int n=mAtoms[o][p];
			xb[n]=xa0[p]-xcm[o];
		}

	}
	for(size_t o=0;o<mCluster.size();o++){
		vector<Dvect> xcm0(mCluster[o].size());
		for(size_t p=0;p<mCluster[o].size();p++){
			xcm0[p]=xcm[mCluster[o][p]];
		}

// Shuffles residue geometric center to get a cluster of minimum size

		auto Encmp=EnComp<T>(xcm0,co);
		cmSweep(xcm0,mCluster[o].size(),Encmp);

// rewrite residue geometric center relative to cluster geometric center

		for(size_t p=0;p<mCluster[o].size();p++){
			xcmC[o]+=xcm0[p];
		}
		xcmC[o]/=static_cast<double>(mCluster[o].size());
		for(size_t p=0;p<mCluster[o].size();p++){
			xcm[mCluster[o][p]]=xcm0[p]-xcmC[o];
		}

	}

// Shuffle cluster geometry center to find the minimal size of the ensemble
	vector<Dvect> xcmC0(mCluster.size());
	for(size_t p=0;p<mCluster.size();p++){
		xcmC0[p]=xcmC[p];
	}
	auto Encmp=EnComp<T>(xcmC0,co);
	cmSweep(xcmC,mCluster.size(),Encmp);
	for(size_t p=0;p<mCluster.size();p++){
		xcmCell+=xcmC0[p];
	}
	xcmCell/=(T) mCluster.size();

// Rewrite the cluster geometric center relative to the aggregate geometric center

	for(size_t p=0;p<mCluster.size();p++){
		xcmC[p]=xcmC0[p]-xcmCell;
	}
	static int UI=0;

// Finally reconstruct the solute reduced atomic positions

	for(size_t o=0;o<mCluster.size();o++){
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			Dvect tmp=xcmCell+xcmC[o]+xcm[n];
			for(size_t i=0;i<mAtoms[n].size();i++){
				int ia=mAtoms[n][i];
				xa[ia]=tmp+xb[ia];
			}
		}

	}

// Place the aggregate in the center of the cell; thus xcmCell is equal to HALF
	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			 xa[ia][o]-=xcmCell[o]-HALF;
			 if(atSolv[ia]) {  // for solvent atoms apply pbc and place to obtain their coordinates in the (0,0,0) cell
				 xa[ia][o]-=rint(xa[ia][o]-HALF);
			 }
		}
	}

	//> Obtain new Cartesian coordinates from reduced xa's
	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			x[ia][o]=Mt.getCO()[o][XX]*xa[ia][XX]+Mt.getCO()[o][YY]*xa[ia][YY]+Mt.getCO()[o][ZZ]*xa[ia][ZZ];
		}
	}
}

template <typename T>
void Atoms<T>::Reconstruct1(Contacts<T> * con0){

	vector<Dvect> nxyz=PBCvect<T>::getVec();
	vector<vector<int> > mCluster=Perco->getCluster();
	vector<vector<int> > mAtoms=Perco->getAtoms();


	Matrix co=Mt.getCO();
	Matrix oc=Mt.getOC();
	double v[DIM];
	v[XX]=co[XX][XX];v[YY]=co[YY][YY];v[ZZ]=co[ZZ][ZZ];
	double MinCO=*min_element(v,v+3)*0.5;  // Half the smallest axis is the cutoff distance for atom shifts
	Dvect Xref;

	// Reconstruct each one of the molecules
	vector<Dvect> xcm(mAtoms.size(),Dvect{T{0.0}});

	for(size_t o=0;o<mAtoms.size();o++){

		Xref=xa[mAtoms[o][0]];
		for(size_t p=0;p<mAtoms[o].size();p++){
			int n=mAtoms[o][p];
			Dvect x1=Xref-xa[n];
			Dvect xp=xa[n];
			Dvect xc1=co*x1;
			if(xc1.Norm() > MinCO) {
				Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),CellComp<T>(x1,co));
				for(int q=0; q < DIM;q++) xa[n][q]=xa[n][q]-tmp[q];
			}
			Xref=xa[n];
		}
		// Compute molecule center of mass xcm
		for(size_t p=0;p<mAtoms[o].size();p++){
			Dvect xx=xa[mAtoms[o][p]];
			xcm[o]+=xx;
		}
		xcm[o]/=static_cast<double>(mAtoms[o].size());
	}
	vector<bool> atSolv(nr,true);
	for(auto o=0;o<mAtoms.size();o++){
		for(auto p=0;p<mAtoms[o].size();p++){
			atSolv[mAtoms[o][p]]=false;
		}
	}

	//
	vector<Dvect> xcmC(mCluster.size());
	auto cmSweep=[this,mCluster,mAtoms,nxyz](RgComp<T> & Rgcmp,vector<Dvect> & xcm0, int o){

		for(size_t p=0;p<mCluster[o].size();p++){
			size_t p0=mCluster[o][p];
			Dvect Xref=xcm0[p];
			Rgcmp.setVect(Xref);
			Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),Rgcmp);
			for(int q=0; q < DIM;q++) xcm0[p][q]=xcm0[p][q]+tmp[q];
			for(size_t i=0;i<mAtoms[p0].size();i++){
			  int ia=mAtoms[p0][i];
			  for(int q=0; q < DIM;q++) xa[ia][q]=xa[ia][q]+tmp[q];
			}
			Rgcmp.CompRg();
		}
		return Rgcmp.getRg();
	};


	for(size_t o=0;o<mCluster.size();o++){
		if(mCluster[o].size() < 5) continue;
		vector<Dvect> xcm0(mCluster[o].size());
		for(size_t p=0;p<mCluster[o].size();p++){
			xcm0[p]=xcm[mCluster[o][p]];
		}
		(*con0)(xcm0,co,oc);
		con0->Neighbors();

		vector<vector<int>> & nnl=con0->NNL();


// Cycle on all center of mass and find the aggregate with minimum Rg;
		auto Rgcmp=RgComp<T>(xcm0,co);
		bool notOk{true};
		T Rg{1e10};
		int M{0};
		while(notOk && M <4){
			T Rg_nw=cmSweep(Rgcmp,xcm0,o);
			if(Rg_nw == Rg){
				notOk=false;
			}
			Rg=Rg_nw;
			M++;
		}

		for(size_t p=0;p<mCluster[o].size();p++){
			xcm[p]=xcm0[p];
		}
		// Compute the center of mass of the cluster
		xcmC[o]=0.0;
		for(size_t p=0;p<mCluster[o].size();p++){
			xcmC[o]+=xcm[p];
		}
		xcmC[o]/=static_cast<double>(mCluster[o].size());
	}
	if(mCluster.size() == 1) return __ReconstructOneCluster(atSolv);

	FindCell<T> myCell(Mt.getCO());
	//> Translate the cell to get all images inside the cell.
	Dvect Translation{myCell.Run(mCluster,mAtoms,xa,nr)};

	//> Translate all coordinates as the center of the super cluster must be the center of the system

	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			xa[ia][o]+=Translation[o];
			if(atSolv[ia]) xa[ia][o]-=rint(xa[ia][o]-HALF);
		}
	}
	for(size_t o=0;o<mAtoms.size();o++){
		// Compute molecule center of mass xcm
		xcm[o]=0.0;
		for(size_t p=0;p<mAtoms[o].size();p++){
			Dvect xx=xa[mAtoms[o][p]];
			xcm[o]+=xx;
		}
		xcm[o]/=static_cast<double>(mAtoms[o].size());
	}
	for(size_t o=0;o<mCluster.size();o++){
		Dvect xcmC{T{0.0}};
		for(size_t p=0;p<mCluster[o].size();p++){
			xcmC+=xcm[mCluster[o][p]];
		}
		xcmC/=static_cast<double>(mCluster[o].size());
		for(size_t p=0;p<mCluster[o].size();p++){
			int n=mCluster[o][p];
			for(size_t i=0;i<mAtoms[n].size();i++){
				int ia=mAtoms[n][i];
				xa[ia][XX]=xa[ia][XX]-rint(xcmC[XX]-HALF);
				xa[ia][YY]=xa[ia][YY]-rint(xcmC[YY]-HALF);
				xa[ia][ZZ]=xa[ia][ZZ]-rint(xcmC[ZZ]-HALF);
			}

		}
	}


	//> Obtain new Cartesian coordinates from reduced xa's
	for(auto ia=0;ia<nr;ia++){
		for(int o=0;o<DIM;o++){
			x[ia][o]=Mt.getCO()[o][XX]*xa[ia][XX]+Mt.getCO()[o][YY]*xa[ia][YY]+Mt.getCO()[o][ZZ]*xa[ia][ZZ];
		}
	}
}

//template <typename T>
//void Atoms<T>::Reconstruct(Contacts<T> * con0){
//
//	vector<Dvect> nxyz=PBCvect<T>::getVec();
//	vector<vector<int> > mCluster=Perco->getCluster();
//	vector<vector<int> > mAtoms=Perco->getAtoms();
//
//	Matrix co=Mt.getCO();
//	Matrix oc=Mt.getOC();
//	T v[DIM];
//	v[XX]=co[XX][XX];v[YY]=co[YY][YY];v[ZZ]=co[ZZ][ZZ];
//
//	T MinCO=*min_element(v,v+3)*0.5;  // Half the smallest axis is the cutoff distance for atom shifts
//	Dvect Xref;
//
//	// Reconstruct each one of the molecules
//	vector<Dvect> xcm(mAtoms.size(),0);
//
//	for(size_t o=0;o<mAtoms.size();o++){
//
//		Xref=xa[mAtoms[o][0]];
//		for(size_t p=0;p<mAtoms[o].size();p++){
//			int n=mAtoms[o][p];
//			Dvect x1=Xref-xa[n];
//			Dvect xp=xa[n];
//			Dvect xc1=co*x1;
//			if(xc1.Norm() > MinCO) {
//				Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),CellComp<T>(x1,co));
//				for(int q=0; q < DIM;q++) xa[n][q]=xa[n][q]-tmp[q];
//			}
//			Xref=xa[n];
//		}
//		// Compute molecule center of mass xcm
//		for(size_t p=0;p<mAtoms[o].size();p++){
//			Dvect xx=xa[mAtoms[o][p]];
//			xcm[o]+=xx;
//		}
//		xcm[o]/=static_cast<T>(mAtoms[o].size());
//	}
//
//	//
//	for(size_t o=0;o<mCluster.size();o++){
//		Dvect xcmC{0};
//		vector<Dvect> xcm0(mCluster[o].size());
//
//		for(size_t p=0;p<mCluster[o].size();p++){
//			xcm0[p]=xcm[mCluster[o][p]];
//		}
//
////		Contacts * con0=new Contacts(xcm0,co,oc);
//
//		(*con0)(xcm0,co,oc);
//		con0->Neighbors();
//		size_t p=0,q=0,qq=0;
//// Start with a reference molecule p=0;
//		Xref=xcm0[p]; // Initial reference molecule is the first on the list
//		p=con0->next(); //p is now the next molecule closest to p-1
//
//		while (p < SIZE_T){
//			size_t n=mCluster[o][p];
//			Dvect x1=Xref-xcm0[p];
//			Dvect xc1=co*x1;
//
//			if(xc1.Norm() >= MinCO) {
//				Dvect tmp=*min_element(nxyz.begin(),nxyz.end(),CellComp<T>(x1,co));
//				for(int q=0; q < DIM;q++) xcm[n][q]=xcm[n][q]-tmp[q];
//				for(int q=0; q < DIM;q++) xcm0[p][q]=xcm0[p][q]-tmp[q];
//
//				for(size_t i=0;i<mAtoms[n].size();i++){
//					int ia=mAtoms[n][i];
//					for(int q=0; q < DIM;q++) xa[ia][q]=xa[ia][q]-tmp[q];
//					}
//
//			}
//
//			Xref=xcm0[p];
//			p=con0->next(); //p is now the next molecule closest to p-1
//
//		};
//		// Compute the center of mass of the cluster
//		for(size_t p=0;p<mCluster[o].size();p++){
//			xcmC+=xcm0[p];
//		}
//		xcmC/=static_cast<T>(mCluster[o].size());
//		for(size_t p=0;p<mCluster[o].size();p++){
//			int n=mCluster[o][p];
//			for(size_t i=0;i<mAtoms[n].size();i++){
//				int ia=mAtoms[n][i];
//				xa[ia][XX]=xa[ia][XX]-rint(xcmC[XX]-HALF);
//				xa[ia][YY]=xa[ia][YY]-rint(xcmC[YY]-HALF);
//				xa[ia][ZZ]=xa[ia][ZZ]-rint(xcmC[ZZ]-HALF);
//			}
//
//		}
//
//	}
//
//	//> Obtain new Cartesian coordinates from reduced xa's
//
//	for(size_t i1=0;i1<mCluster.size();i1++){
//		for(size_t i2=0;i2<mCluster[i1].size();i2++){
//			int n=mCluster[i1][i2];
//			for(size_t i=0;i<mAtoms[n].size();i++){
//				int ia=mAtoms[n][i];
//				for(int o=0;o<DIM;o++){
//					xa[ia][o]=xa[ia][o]+0.5;
//					x[ia][o]=Mt.getCO()[o][XX]*xa[ia][XX]+Mt.getCO()[o][YY]*xa[ia][YY]+Mt.getCO()[o][ZZ]*xa[ia][ZZ];
//				}
//			}
//		}
//	}
//}
template <typename T>
DDvect<T> Atoms<T>::__FindCell(const vector<vector<int> > & mCluster, const vector<vector<int> > & mAtoms){
	vector<Dvect> x(nr);
	Matrix CO(Mt.getCO());
	for(auto o=0;o<nr;o++){
		for(auto n=0;n<DIM;n++){
			x[o][n]=xa[o][n];
		}
	}
	vector<T> xcm(mAtoms.size(),0.0);
	Dvect myReturn{T{0.0}};
	T dx=0.2;
	for(auto n=0;n<DIM;n++){
		size_t Nn=CO[n][n]/dx;
		size_t M{0};
		T ddx=1.0/Nn;
		vector<int> myPBC(Nn);
		while(M < Nn){
			size_t bCount{0};
			for(auto o=0;o<nr;o++){
				x[o][n]+=ddx;
			}
			for(size_t o=0;o<mAtoms.size();o++){
				// Compute molecule center of mass xcm
				xcm[o]=0.0;
				for(size_t p=0;p<mAtoms[o].size();p++){
					T xx=x[mAtoms[o][p]][n];
					xcm[o]+=xx;
				}
				xcm[o]/=static_cast<T>(mAtoms[o].size());
				xcm[o]-=rint(xcm[o]-0.5);
			}

			for(size_t o=0;o<mCluster.size();o++){
				int nn=mCluster[o][0];
				T Xref=xcm[nn];
				for(size_t p=1;p<mCluster[o].size();p++){
					int nn=mCluster[o][p];
					T PBC=rint(Xref-xcm[nn]);
					if(PBC) bCount++;
				}
			}
			myPBC[M]=bCount;
			M++;
		}
		auto result = std::minmax_element(myPBC.begin(),myPBC.end());
		auto it1=std::find_if(myPBC.begin(),myPBC.end(),[result](int q){return q==*result.first;});
		auto it2=std::find_if(it1,myPBC.end(),[result](int q){return q!=*result.first;});
		auto j=std::distance(myPBC.begin(),it1)+std::distance(it1,it2)/2;
		myReturn[n]=j*ddx;
	}
	return 	myReturn;
}
template <typename T>
Atoms<T>::~Atoms(){
	delete Perco;
	for(auto & op: Rg_i){
		delete op;
	}

}
template class Atoms<float>;
template class Atoms<double>;
template class CellComp<float>;
template class CellComp<double>;
template class RgComp<float>;
template class RgComp<double>;
template void Atoms<double>::Gyro<Enums::noJSON>();
template void Atoms<double>::Gyro<Enums::JSON>();
template void Atoms<float>::Gyro<Enums::noJSON>();
template void Atoms<float>::Gyro<Enums::JSON>();
template void Atoms<double>::SetupPercolate<Enums::noJSON>(Topol_NS::Topol &);
template void Atoms<double>::SetupPercolate<Enums::JSON>(Topol_NS::Topol &);
template void Atoms<float>::SetupPercolate<Enums::noJSON>(Topol_NS::Topol &);
template void Atoms<float>::SetupPercolate<Enums::JSON>(Topol_NS::Topol &);
template void Atoms<float>::InitSelection<Enums::Selection>(vector<string> &, Topol_NS::Topol &);
template void Atoms<float>::InitSelection<Enums::Reference>(vector<string> &, Topol_NS::Topol &);
template void Atoms<float>::InitSelection<Enums::fftPadding>(vector<string> &, Topol_NS::Topol &);
template void Atoms<double>::InitSelection<Enums::Selection>(vector<string> &, Topol_NS::Topol &);
template void Atoms<double>::InitSelection<Enums::Reference>(vector<string> &, Topol_NS::Topol &);
template void Atoms<double>::InitSelection<Enums::fftPadding>(vector<string> &, Topol_NS::Topol &);

