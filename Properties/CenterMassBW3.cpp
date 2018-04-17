/*
 * CenterMassBW3.cpp
 *
 *  Created on: May 12, 2015
 *      Author: marchi
 */

#include "CenterMassBW3.h"

template <typename T>
bool CenterMassBW3<T>::BeenThere=false;

template <typename T>
void CenterMassBW3<T>::SkipIt(ifstream & fin){
	if(AlreadySkipped){
		fin.seekg(SkipPos,ios::cur);
		return;
	}
	int r;
	size_t mm,nn;
	T t;
	size_t nQm;
	streampos pos0,pos1;
	pos0=fin.tellg();
	fin.read(reinterpret_cast<char *> (&time_c), sizeof(time_c));
	fin.read(reinterpret_cast<char *> (&r), sizeof(r));
	vector<int> rr(r);
	fin.read(reinterpret_cast<char *> (&rr[0]), (sizeof rr[0])*r);
	fin.read(reinterpret_cast<char *> (&nn), sizeof(nn));
	fin.read(reinterpret_cast<char *> (&mm), sizeof(mm));

	streampos pos=static_cast<streampos> (2*sizeof(Matrix)+sizeof(t)*(r*DIM+r*(DIM+1)+2*nn*DIM));
	fin.seekg(pos,ios::cur);
	if(!Labels.size()) {
		Labels=vector<string>(mm);
		for(auto p=0;p<mm;p++) {
			string labels(5,' ');
			fin.read(&labels[0],sizeof(char)*5);
			labels.erase(remove_if(labels.begin(),labels.end(),::isspace),labels.end());
			Labels[p]=labels;
		}
	}
	pos=static_cast<streampos> (sizeof(t)*(2*mm*DIM));
	fin.seekg(pos,ios::cur);
	fin.read(reinterpret_cast<char *> (&nQm), (sizeof nQm));
	pos=static_cast<streampos>(sizeof(t)*(nQm*(DIM+1)));
	fin.seekg(pos,ios::cur);
	pos1=fin.tellg();
	SkipPos=pos1-pos0-static_cast<streampos>(5*mm*sizeof(char));
	AlreadySkipped=true;
}
template <typename T>
void CenterMassBW3<T>::ReadIt(ifstream & fin){
	auto & Dt=this->Dt;
	auto & BeenThere=this->BeenThere;
	auto & CO=this->CO;
	auto & COt=this->COt;
	auto & OC=this->OC;
	auto & OCt=this->OCt;
	auto & cmt=this->cmt;
	auto & Q_tt=this->Q_tt;
	auto & Qmt=this->Qmt;

	int r;
	size_t mm,nn;
	CenterMass_t PickDiff=CMPick<T>()(this);

	fin.read(reinterpret_cast<char *> (&time_c), sizeof(time_c));

	if(Dt < 0.0) Dt+=time_c;
	if(!BeenThere) Dt=-time_c;

	fin.read(reinterpret_cast<char *> (&r), sizeof(r));
	vector<int> rr(r);
	fin.read(reinterpret_cast<char *> (&rr[0]), (sizeof rr[0])*r);
	fin.read(reinterpret_cast<char *> (&nn), sizeof(nn));
	fin.read(reinterpret_cast<char *> (&mm), sizeof(mm));
	fin.read(reinterpret_cast<char *> (&CO), sizeof(Matrix));
	fin.read(reinterpret_cast<char *> (&OC), sizeof(Matrix));
	if(PickDiff == diffMicelles){
		vector<T> t(r*DIM);
		vector<T> v(r*(DIM+1));
		fin.read(reinterpret_cast<char *> (&t[0]), sizeof(t[0])*r*DIM);
		fin.read(reinterpret_cast<char *> (&v[0]), sizeof(v[0])*r*(DIM+1));
		if(!BeenThere) {
			cmt=vector<vector<Dvect> >(r);
			Q_tt=vector<vector<Quaternion> >(r);
		}
		for(auto o=0;o<r;o++){
			Dvect myt{t[DIM*o],t[DIM*o+1],t[DIM*o+2]};
			Quaternion myv{v[(DIM+1)*o],v[(DIM+1)*o+1],v[(DIM+1)*o+2],v[(DIM+1)*o+3]};
			cmt[o].push_back(myt);
			Q_tt[o].push_back(myv);
		}

	}else{
		streampos pos=static_cast<streampos>(sizeof(T)*(r*DIM+r*(DIM+1)));
		fin.seekg(pos,ios::cur);
	}
	if(true){
		streampos pos=static_cast<streampos>(sizeof(T)*(2*nn*DIM));
		fin.seekg(pos,ios::cur);
	} else{
		vector<T> tmR(nn*DIM);
		vector<T> tmL(nn*DIM);
		fin.read(reinterpret_cast<char *>(&tmR[0]),sizeof(tmR[0])*nn*DIM);
		fin.read(reinterpret_cast<char *>(&tmL[0]),sizeof(tmL[0])*nn*DIM);
	}

	if(!Labels.size()) {
		Labels=vector<string>(mm);
		for(auto p=0;p<mm;p++) {
			string labels(5,' ');
			fin.read(&labels[0],sizeof(char)*5);
			labels.erase(remove_if(labels.begin(),labels.end(),::isspace),labels.end());
			Labels[p]=labels;
		}
	}

	if(PickDiff == diffk){
		if(!BeenThere) {
			X=vvector_d(mm,vector<Dvect>());
			Xd=vvector_d(mm,vector<Dvect>());
		}
		vector<T> tmX(mm*DIM);
		vector<T> tmXd(mm*DIM);
		fin.read(reinterpret_cast<char *> (&tmX[0]), (sizeof tmX[0])*mm*DIM);   // In Cartesian coordinates
		fin.read(reinterpret_cast<char *> (&tmXd[0]), (sizeof tmXd[0])*mm*DIM); // In Cartesian coordinates
		vector<Dvect> x(mm);

		for(auto n=0;n<mm;n++){
			for(auto o=0;o<DIM;o++)
				x[n][o]=tmX[DIM*n+o];
			X[n].push_back(x[n]);
			for(auto o=0;o<DIM;o++)
				x[n][o]=tmXd[DIM*n+o];
			Xd[n].push_back(x[n]);
		}
	} else{
		streampos pos=static_cast<streampos>(sizeof(T)*(2*mm*DIM));
		fin.seekg(pos,ios::cur);
	}

	COt.push_back(CO);
	OCt.push_back(OC);
	size_t nQm;
	fin.read(reinterpret_cast<char *> (&nQm), (sizeof nQm));

	if(PickDiff == diffq){
		vector<T> vh((DIM+1)*nQm);
		fin.read(reinterpret_cast<char *> (&vh[0]), (sizeof vh[0])*nQm*(DIM+1));
		if(!BeenThere)
			Qmt=vector<vector<Quaternion> >(nQm);

		for(auto o=0;o<nQm;o++)
			Qmt[o].push_back({vh[(DIM+1)*o],vh[(DIM+1)*o+1],vh[(DIM+1)*o+2],vh[(DIM+1)*o+3]});
	} else{
		streampos pos=static_cast<streampos>(sizeof(T)*(nQm*(DIM+1)));
		fin.seekg(pos,ios::cur);
	}

	BeenThere=true;
}
template <typename T>
vector<complex<T>> CenterMassBW3<T>::DiffQQ(vvector_q & Q,const int ell, const int mp, const int mm){
	size_t Ndim=Q[0].size();
	vector<complex<T> > Ctot(Ndim,{0.0,0.0});
	for(auto pv: Q){
		alglib::complex_1d_array U,Uc,C;
		U.setlength(Ndim);
		C.setlength(2*Ndim-1);
		Uc.setlength(Ndim);
		for(auto o=0;o<Ndim;o++){
			WignerDMatrix W;
			W.SetRotation(pv[o]);
			auto s=W(ell,mp,mm);
			alglib::complex t(s.real(),s.imag());
			U[o]={s.real(),s.imag()};
			Uc[o]=U[o];
		}
		alglib::corrc1d(U,Ndim,Uc,Ndim,C);
		for(auto o=0;o<Ndim;o++){
			C[o].x/=static_cast<T> (Ndim-o);
			C[o].y/=static_cast<T> (Ndim-o);
		}

		for(auto o=0;o<Ndim;o++){
			T x=static_cast<T>(C[o].x);
			T y=static_cast<T>(C[o].y);
			Ctot[o]+=complex<T>{x,y};
		}
	}

	for(auto o=0;o<Ndim;o++){
		Ctot[o]/=static_cast<T> (Q.size());
	}
	return Ctot;
}


template <typename T>
alglib::real_1d_array CenterMassBW3<T>::Diff1D(vector<Dvect> & X, size_t Ndim){
	alglib::real_1d_array Ux,Uy,Uz,Cx,Cy,Cz,Sab,Saa,S;
	Ux.setlength(Ndim);
	Uy.setlength(Ndim);
	Uz.setlength(Ndim);
	Cx.setlength(2*Ndim-1);
	Cy.setlength(2*Ndim-1);
	Cz.setlength(2*Ndim-1);
	Sab.setlength(Ndim);
	Saa.setlength(Ndim);
	S.setlength(Ndim);
	for(auto o=0;o<X.size();o++){
		Ux[o]=X[o][XX];
		Uy[o]=X[o][YY];
		Uz[o]=X[o][ZZ];
	}
	T dsq=0.0;
	for(auto o=0;o<Ndim;o++)
		dsq+=Ux[o]*Ux[o]+Uy[o]*Uy[o]+Uz[o]*Uz[o];


	alglib::corrr1d(Ux,Ndim,Ux,Ndim,Cx);
	alglib::corrr1d(Uy,Ndim,Uy,Ndim,Cy);
	alglib::corrr1d(Uz,Ndim,Uz,Ndim,Cz);
	for(auto o=0;o<Ndim;o++) Sab[o]=Cx[o]+Cy[o]+Cz[o];
	Saa[0]=2.0*dsq;
	for(auto o=1;o<Ndim;o++){
		T r0=Ux[o-1]*Ux[o-1]+Uy[o-1]*Uy[o-1]+Uz[o-1]*Uz[o-1];
		T r1=Ux[Ndim-o]*Ux[Ndim-o]+Uy[Ndim-o]*Uy[Ndim-o]+Uz[Ndim-o]*Uz[Ndim-o];
		Saa[o]=Saa[o-1]-r0-r1;
	}
	for(auto o=0;o<Ndim;o++){
		S[o]=(Saa[o]-2.0*Sab[o])/static_cast<T>(Ndim-o);
	}
	return S;
}
template <typename T>
void CenterMassBW3<T>::Diffusion(ofstream & fout){
	auto & Dt=this->Dt;

	alglib::real_1d_array S,Sd;
	size_t Ndim=X[0].size();

	S.setlength(Ndim);
	Sd.setlength(Ndim);
	for(auto it=Labels.begin();it != Labels.end();++it){
		if(Diffs.find(*it) == Diffs.end()){
			Diffs[*it]=vector<dual_d<T>>(Ndim,{0.0,0.0});
			nDiffs[*it]=0;
		}
	}
	for(auto n=0;n<X.size();n++){
		S=Diff1D(X[n],Ndim);
		Sd=Diff1D(Xd[n],Ndim);
		for(auto o=0;o<Ndim;o++) {
			Diffs[Labels[n]][o].x+=S[o];
			Diffs[Labels[n]][o].y+=Sd[o];
		}
		nDiffs[Labels[n]]++;

	}
	auto nLegend=0;
	fout << "# Grace file " <<endl;
	fout << "# Times are in ns" <<endl;
	fout << "# Distances are in nm^2" <<endl;
	fout << "# " <<endl;
	for(auto it=Diffs.begin();it != Diffs.end();++it){
		long int accu=nDiffs[it->first];
		for(auto o=0;o<Ndim;o++){
			it->second[o].x/=static_cast<T>(accu);
			it->second[o].y/=static_cast<T>(accu);
		}
		stringstream legend;
		legend << "@    S";
		legend.width(3);
		legend << std::fixed << std::left << nLegend++ ;
		legend<< " legend   "<< "\"" << it->first << "\"" <<endl;
		fout << legend.str() ;
	}
	stringstream legend;
	legend << "@    S";
	legend.width(3);
	legend << std::fixed << std::left << nLegend++ ;
	legend<< " legend   "<< "\"" << " Rot " << "\"" <<endl;
	fout << legend.str() ;
	nLegend=0;
	for(auto it=Diffs.begin();it != Diffs.end();++it){
		stringstream legend;
		legend << "@target G0.S";
		legend.width(3);
		legend << std::fixed << std::left << nLegend++ << endl;
		legend << "@type xy" <<endl;
		fout << legend.str();
		for(auto o=0;o< it->second.size();o++){
			fout << std::fixed<< std::setprecision(3) << Dt*static_cast<T>(o)/1000.0 ;
			fout << " " ;
			fout << std::fixed << std::setw(10)<< std::setprecision(5)<< it->second[o].x << endl;
		}
		fout << "&"<<endl;
	}
	legend.clear();
	legend << "@target G0.S";
	legend.width(3);
	legend << std::fixed << std::left << nLegend++ << endl;
	legend << "@type xy" <<endl;
	fout << legend.str();

	for(auto o=0;o< Diffs.begin()->second.size();o++){
		fout << std::fixed<< std::setprecision(3) << Dt*static_cast<T>(o)/1000.0 ;
		fout << " " ;
		fout << std::fixed << std::setw(10)<< std::setprecision(5)<< Diffs.begin()->second[o].y << endl;
	}
	fout << "&"<<endl;

}
template <typename T>
void CenterMassBW3<T>::WriteIt(ostream & fout) {
	auto & CO=this->CO;
	auto & OC=this->OC;
	auto & cm=this->cm;
	auto & cmm=this->cmm;
	auto & Q_t=this->Q_t;
	vector<T> t(cm.size()*DIM);
	vector<T> v;
	for(auto q: Q_t){
		Quaternion qq=q.normalized();
		v.push_back(qq[0]);
		v.push_back(qq[1]);
		v.push_back(qq[2]);
		v.push_back(qq[3]);
	}
	auto n=0;
	for(auto g0: cm)
		for(auto q=0;q<DIM;q++)
			t[n++]=g0[q];

	int r=cm.size();
	vector<int> rr(r);
	for(auto i=0;i<r;i++) rr[i]=cmm[i].size();

	vector<T> tmR;
	vector<T> tmL;
	vector<T> tmX;
	vector<T> tmXd;


	double time=CenterMass<T>::getTime();
	vector<vector<Dvect> > cmL=cmm;
	vector<vector<Dvect> > cmR=cmm;
	this->RotateBack(cmR);

	auto po=[](vector<T> & x,vector<vector<Dvect> > & y){
		auto p=0;for(auto & R0: y)
			for(auto & R1: R0)
				for(auto m=0;m<DIM;m++)
					x.push_back(R1[m]); };
	po(tmR,cmR);
	po(tmL,cmL);
	po(tmX,xc);
	po(tmXd,xd);
	auto nn=tmR.size()/DIM;
	auto mm=tmX.size()/DIM;
	fout.write((char *) &time, sizeof time);
	fout.write((char *) &r, sizeof r);
	fout.write((char *) &rr[0], (sizeof rr[0])*r);
	fout.write((char *) &nn, sizeof nn);
	fout.write((char *) &mm, sizeof mm);
	fout.write(reinterpret_cast<const char *> (&CO), sizeof(Matrix));
	fout.write(reinterpret_cast<const char *> (&OC), sizeof(Matrix));
	fout.write((char *) &t[0], (sizeof t[0])*r*DIM);
	fout.write((char *) &v[0], (sizeof v[0])*r*(DIM+1));
	fout.write((char *) &tmR[0], (sizeof tmR[0])*nn*DIM);
	fout.write((char *) &tmL[0], (sizeof tmL[0])*nn*DIM);

	if(!BeenThere) {
		const size_t strDIM=5;
		for(auto o=0;o<mm;o++){
			Labels[o].append(string(strDIM,' '),Labels[o].size()-1,strDIM-Labels[o].size());
			fout.write(Labels[o].c_str(),sizeof(char)*strDIM);
		}
	}
	fout.write(reinterpret_cast<const char *> (&tmX[0]), (sizeof tmX[0])*mm*DIM);
	fout.write(reinterpret_cast<const char *> (&tmXd[0]), (sizeof tmXd[0])*mm*DIM);
	size_t nQm=this->Qm.size();
	v.clear();
	for(auto q: this->Qm){
		Quaternion qq=q.normalized();
		v.push_back(qq[0]);
		v.push_back(qq[1]);
		v.push_back(qq[2]);
		v.push_back(qq[3]);
	}
	fout.write(reinterpret_cast<char *> (&nQm), (sizeof nQm));
	fout.write(reinterpret_cast<char *> (&v[0]), (sizeof v[0])*nQm*(DIM+1));

	BeenThere=true;
}


template <typename T>
CenterMassBW3<T>::CenterMassBW3() {
}
template <typename T>
CenterMassBW3<T>::~CenterMassBW3() {
	// TODO Auto-generated destructor stub
}


template <typename T>
string CenterMassBW3<T>::HeavyAtoms1="S C";
template <typename T>
string CenterMassBW3<T>::HeavyAtoms2="O S C";

template <typename T>
void CenterMassBW3<T>::RotateTranslate(vvector_d & t, vvector_d & s){
	s=vvector_d(t.size());
	for(auto n=0;n< t.size();n++){
		auto qq=this->Q_t[n];
		auto Rot=this->MatrixfromQ(qq.inverse());
		auto Rot1=this->MatrixfromQ(qq);
		auto R_cm=this->OC*this->cm[n];
		for(auto & pv: t[n]){
			Dvect t0=pv;
			Dvect t1=t0-R_cm;
			for(auto o=0;o<DIM;o++) t1[o]=t1[o]-rint(t1[o]);
			t0=this->CO*t1;
			Dvect yy={R00,0.0,0.0};
			pv=t0;
			t0=Rot1*yy;
			s[n].push_back(t0);
		}
	}
}

template <typename T>
void CenterMassBW3<T>::RefCoord(rvec * xa, vector<string> & atn, vvector_i & ind, vvector_i & mats){
	if(Labels.empty()) {
		for(auto o=0;o<ind.size();o++)
			for(auto p=0;p<ind[o].size();p++)
				for(auto q=0;q<mats[p].size();q++){
					auto r=mats[p][q];
					if(HeavyAtoms1.find(atn[r][0]) != string::npos)
						Labels.push_back(atn[r]);
				}
	}
	vvector_d t(ind.size());
	vvector_d s(ind.size());
	for(auto o=0;o<ind.size();o++){
		for(auto p=0;p<ind[o].size();p++){
			auto op=ind[o][p];
			for(auto q=0;q<mats[op].size();q++){
				auto r=mats[op][q];
				if(HeavyAtoms1.find(atn[r][0]) != string::npos){
					Dvect tu;
					tu[XX]=static_cast<T>(xa[r][XX]);
					tu[YY]=static_cast<T>(xa[r][YY]);
					tu[ZZ]=static_cast<T>(xa[r][ZZ]);
					t[o].push_back(tu);
				}

			}
		}
	}

	RotateTranslate(t,s);
	xd=s;
	xc=t;

}

template <typename T>
void CenterMassBW3<T>::operator ()(vector<Dvect> & v,  vector<vector<Dvect> > & vm, vector<Quaternion> & Q, vector<Quaternion> & Q_m) {
	CenterMass<T>::operator()(v,Q);
	for(auto i=0;i<vm.size();i++){
		try{if(this->DoneOnce) if(this->cmm[i].size() != vm[i].size()) throw "The number of molecules in clusters changes from step to step!! ";}
		catch(const char * s){cout << s <<endl;exit(1);}
		this->cmm[i]=vm[i];
	}

	try{if(this->DoneOnce) if(Q_m.size() != this->Qm.size()) throw "The number of molecular quaternions in clusters changes from step to step!! ";}
	catch(const char * s){cout << s <<endl;exit(1);}

	for(auto n=0;n<Q_m.size();n++)
		this->Qm[n]=this->Qm[n]*Q_m[n];
	this->DoneOnce=true;
}
template class CenterMassBW3<float>;
template class CenterMassBW3<double>;

