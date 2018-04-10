/*
 * CenterMassWigner.cpp
 *
 *  Created on: Jun 1, 2015
 *      Author: marchi
 */

#include "CenterMassWigner.h"
//enum diffAction_t {diffk, diffq, diffMicelles};

template <typename T>
void CenterMassWigner<T>::DiffQ(ofstream & fout,const int ell, const int mp, const int mm){

	size_t Ndim=this->Qmt[0].size();
	vector<complex<T> > Ctot=this->DiffQQ(this->Qmt,ell,mp,mm);

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
		fout << std::fixed<< std::setprecision(3) << this->Dt*static_cast<double>(o)/1000.0 ;
		fout << " " ;
		fout << std::fixed << std::setw(10)<< std::setprecision(5)<< Ctot[o].real()/Ctot[0].real() << " ";
		fout << std::fixed << std::setw(10)<< std::setprecision(5)<< Ctot[o].imag() <<endl; ;
	}

}

template <typename T>
CenterMassWigner<T>::CenterMassWigner() {
	// TODO Auto-generated constructor stub

}

template <typename T>
CenterMassWigner<T>::~CenterMassWigner() {
	// TODO Auto-generated destructor stub
}
template class CenterMassWigner<float>;
template class CenterMassWigner<double>;

