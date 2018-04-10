/*
 * VoronoiBinary.cpp
 *
 *  Created on: Dec 2, 2017
 *      Author: marchi
 */

#include "VoronoiBinary.h"

namespace Voro {
template <typename T>
void VoronoiBinary<T>::WriteIt(std::ofstream & fout){
	if(this->writeBinary)
		this->bPrintBody(fout);
	else
		T::WriteIt(fout);

}
template <typename T>
void VoronoiBinary<T>::ReadIt(std::ifstream & fin){
	this->bReadBody(fin);
}

template <typename T>
VoronoiBinary<T>::VoronoiBinary(ifstream & fin) {
	this->bReadHeader(fin);
	this->Vol=vector<double>(this->nr);
	this->Neighs=vector<vector<int>>(this->nr);
	this->Surface=vector<vector<double>>(this->nr);
	this->area.Allocate(this->nresid,this->nc);
	this->Vols.Allocate(this->nresid);
	this->wShells=vector<vector<int>>(VoronoiSetter::maxLevel);

	this->readBinary=true;
}
template <typename T>
VoronoiBinary<T>::VoronoiBinary(ofstream & fout,Topol & myTop,bool bH, Parallel::NewMPI * curr): T(myTop,bH){
	if(!curr->Get_Rank())
		this->bPrintHeader(fout);
	this->writeBinary=true;
}

template <typename T>
VoronoiBinary<T>::~VoronoiBinary() {
	// TODO Auto-generated destructor stub
}
template class VoronoiBinary<VoronoiMicelles>;
template class VoronoiBinary<VoronoiMicellesJSON>;
} /* namespace Voro */
