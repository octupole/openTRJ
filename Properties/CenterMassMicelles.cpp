/*
 * CenterMassMicelles.cpp
 *
 *  Created on: Jun 3, 2015
 *      Author: marchi
 */

#include "CenterMassMicelles.h"

Histogram1Db::Histogram1Db(double dx0, double cut0): Histogram1D() {
		cut=cut0;
		dx=dx0;
		Label=" ";
		HisX=static_cast<size_t>(cut/dx);
		hist=vector<hist1D>(HisX+1);
}
void Histogram1Db::Normalize(double rho=1){
		double Norm=0.0;
		for(auto o=0;o<HisX;o++){
			double r=dx*(o+1);
			double Voldr=4.0*PI*r*r*dx;
			hist[o].hisX=(hist[o].gHisx()/Voldr)/rho;
		}
	}
size_t Histogram1Db::size(){
	return hist.size();
}
hist1D & Histogram1Db::operator[](size_t n){
	try{
		if(n > hist.size()) throw "Something wrong with histogram ";
	} catch(const char * s){
		cout << s << " " << n << " " << hist.size() <<endl;
		exit(1);
	}
	return hist.at(n);
}

ostream & operator<<(ostream & fout, Histogram1Db & y){
	fout << "# Radial density  " + y.Label + " "<< endl;
	fout << "#  R       GofR(r) " << endl;
	for(int o=0;o< y.HisX;o++){
		double ddx=y.dx*static_cast<double>(o);
		double f=y[o].gHisx();

		if(!f) continue;
		fout << fixed << setw(8) << setprecision(2) << ddx/unit_nm;
		fout << fixed << setw(12) << right << scientific << setprecision(4) << f;
		fout << endl;
	}
	return fout;
}

template <typename T>
CenterMassMicelles<T>::CenterMassMicelles() {
	// TODO Auto-generated constructor stub

}

template <typename T>
void CenterMassMicelles<T>::DiffQ(ofstream & fout,const int ell, const int mp, const int mm){
	auto& Dt=this->Dt;
	auto& Q_tt=this->Q_tt;
	vector<complex<T> > Ctot=DiffQQ(Q_tt,ell,mp,mm);
	size_t Ndim=Ctot.size();
	auto nLegend=0;
	fout << "# Grace file " <<endl;
	fout << "# Times are in ns" <<endl;
	fout << "# Distances are in nm^2" <<endl;
	fout << "# " <<endl;
	stringstream legend;
	legend << "@    S";
	legend.width(3);
	legend << std::fixed << std::left << nLegend++ ;
	legend<< " legend   "<< "\" Wigner Correlation \"" <<endl;
	fout << legend.str() ;
	for(auto o=0;o<Ndim;o++){
		fout << std::fixed<< std::setprecision(3) << Dt*static_cast<T>(o)/1000.0 ;
		fout << " " ;
		fout << std::fixed << std::setw(10)<< std::setprecision(5)<< Ctot[o].real()/Ctot[0].real() << " ";
		fout << std::fixed << std::setw(10)<< std::setprecision(5)<< Ctot[o].imag() <<endl; ;
	}
	fout << "& " <<endl;
}

template <typename T>
void CenterMassMicelles<T>::CorrectCoords(vvector_d & y){
	auto & OCt=this->OCt;
	auto & COt=this->COt;
	size_t N=y.size();
	size_t Ndim=y[0].size();
	vector<vector<Dvect> > cm(y.size(),vector<Dvect>(y[0].size(),{0.0,0.0,0.0}));
	for(auto o=0;o<N;o++){
		for(auto p=0;p<Ndim;p++){
			cm[o][p]=OCt[p]*y[o][p];
		}
	}
	for(auto o=0;o<N;o++){
		Dvect Dx={0.0,0.0,0.0};
		for(auto p=1;p<Ndim;p++){
			cm[o][p]=cm[o][p]-Dx;
			auto pb=p-1;
			Dvect dd=cm[o][p]-cm[o][pb];
			for(auto q=0;q<DIM;q++)
				dd[q]=rint(dd[q]);
			cm[o][p]=cm[o][p]-dd;
			Dx+=dd;
		}
	}
	for(auto o=0;o<N;o++){
		for(auto p=0;p<Ndim;p++){
			y[o][p]=COt[p]*cm[o][p];
		}
	}
}

template <typename T>
void CenterMassMicelles<T>::GofR(ofstream & fout ){
	auto & OCt=this->OCt;
	auto & COt=this->COt;
	auto & cmt=this->cmt;
	double dx=0.05;
	double rcut=COt[0][XX][XX]*0.5;
	Histogram1Db gofr(dx,rcut);

	gofr.setLabel("GOFR");
	size_t Ntime=cmt[0].size();
	for(auto o=0;o<cmt.size();o++){
		for(auto p=0;p<cmt.size();p++){
			if(p ==o) continue;
			for(auto q=0;q<Ntime;q++){
				Dvect xa=OCt[q]*(cmt[o][q]-cmt[p][q]);
				xa=xa-Rint(xa);
				Dvect xc=COt[q]*xa;
				double dist=xc.Norm();
				if(dist <rcut){
					size_t h=static_cast<size_t>(dist/dx);
					gofr[h]+=hist1D(1.0,0);
				}
			}
		}
	}
	for(auto o=0;o<gofr.size();o++)
		gofr[o]/=static_cast<double>(Ntime*cmt.size());
	double rho=cmt.size()/(COt[0][XX][XX]*COt[0][YY][YY]*COt[0][ZZ][ZZ]);
	gofr.Normalize(rho);
	fout << gofr;
}


template <typename T>
void CenterMassMicelles<T>::Diffusion(ofstream & fout){
	auto & cmt=this->cmt;
	auto & Dt=this->Dt;

	alglib::real_1d_array S;
	size_t Ndim=cmt[0].size();
	size_t N=cmt.size();
	S.setlength(Ndim);
	vector<double> DiffR(Ndim,0.0);
	vector<vector<double> > DiffE(cmt.size(),DiffR);
	CorrectCoords(cmt);
	for(auto o=0;o<Ndim;o++)
		cout << o << " " << cmt[0][o] <<endl;

	for(auto n=0;n<cmt.size();n++){

		S=this->Diff1D(cmt[n],Ndim);
		for(auto o=0;o<Ndim;o++) {
			DiffR[o]+=S[o];
			DiffE[n][o]=S[o];
		}
	}
	for(auto o=0;o<Ndim;o++)
		DiffR[o]/=static_cast<double>(cmt.size());

	auto nLegend=0;
	fout << "# Grace file " <<endl;
	fout << "# Times are in ns" <<endl;
	fout << "# Distances are in nm^2" <<endl;
	fout << "# " <<endl;
	stringstream legend;
	legend << "@    S";
	legend.width(3);
	legend << std::fixed << std::left << nLegend ;
	legend<< " legend   "<< "\"" << "DiffK" << std::fixed << std::left << nLegend++ <<"\"" <<endl;

	fout << legend.str() ;
	for(auto p=0;p<N;p++){
		stringstream legend;
		legend << "@    S";
		legend.width(3);
		legend << std::fixed << std::left << nLegend ;
		legend<< " legend   "<< "\"" << "DiffK" << std::fixed << std::left << nLegend++ << "\"" <<endl;
		fout << legend.str() ;
	}

	for(auto o=0;o<Ndim;o++){
		fout << std::fixed<< std::setprecision(3) << Dt*static_cast<double>(o)/1000.0 ;
		fout << " " ;
		fout << std::fixed << std::setw(10)<< std::setprecision(5)<< DiffR[o] << endl;
	}
	fout << "&"<<endl;
	for(auto p=0;p<N;p++){
		for(auto o=0;o<Ndim;o++){
			fout << std::fixed<< std::setprecision(3) << Dt*static_cast<double>(o)/1000.0 ;
			fout << " " ;
			fout << std::fixed << std::setw(10)<< std::setprecision(5)<< DiffE[p][o] << endl;
		}
		fout << "&"<<endl;
	}
}

template <typename T>
CenterMassMicelles<T>::~CenterMassMicelles() {
	// TODO Auto-generated destructor stub
}

template class CenterMassMicelles<float>;
template class CenterMassMicelles<double>;
