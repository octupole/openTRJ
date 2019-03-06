/*
 * AtomsProp.cpp
 *
 *  Created on: Mar 1, 2019
 *      Author: marchi
 */

#include "AtomsProp.h"

template <typename T>
void AtomsProp<T,radial>::doProperty(){
	vector<vector<int> > mCluster=this->Perco->getCluster();
	vector<vector<int> > mAtoms=this->Perco->getAtoms();
	Matrix co=this->Mt.getCO();
	Matrix oc=this->Mt.getOC();
	Dvect xcmCell{0}; // Geometric center of the aggregate composed of clusters
	vector<Dvect> xcmC(mCluster.size()); // Geometric center of the clusters
	vector<Dvect> xcm(mAtoms.size(),Dvect{T{0.0}}); // Geometric center of the solute residues
	vector<Dvect> xb(this->nr,Dvect{T{0.0}}); // Atomic coordinates relative to the residue geometric center

	vector<bool> atSolv(this->nr,true);

	// find atoms which are solvent

	for(auto o=0;o<mAtoms.size();o++){
		for(auto p=0;p<mAtoms[o].size();p++){
			atSolv[mAtoms[o][p]]=false;
		}
	}

	vector<Gyration<T> *> Rg=vector<Gyration<T>*>(mCluster.size());
	for(auto & ip: Rg)
		ip=new Gyration<T>();
	Gyration<T>::setTime(this->time_c);

	this->CalcGyro(this->mass,Rg);
	for(size_t o{0};o<Rg.size();o++){
		double Rh{5*Rg[o]->gRadg()/3};
		double rcut=rCut+sqrt(Rh);
		xcm[o]=Rg[o]->gXcm();
	}
	for(size_t o{0};o<mCluster.size();o++){
		for(size_t p{0};p<mCluster[o].size();p++){
			auto n=mCluster[o][p];
			for(size_t q{0};q<mAtoms[n].size();q++){
				auto m=mAtoms[n][q];
				auto x_a=xa[m]-xcm[o];


			}
		}

	}
}
template class AtomsProp<float,radial>;
template class AtomsProp<double,radial>;
