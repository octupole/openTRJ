/*
 * AtomsProp.cpp
 *
 *  Created on: Mar 1, 2019
 *      Author: marchi
 */

#include "AtomsProp.h"

template <typename T>
void AtomsProp<T,radial>::Reduce(Parallel::NewMPI * y){
	for(auto it=histograms.begin();it != histograms.end();it++){
		it->second->Reduce(y);
	}
}


template <typename T>
void * AtomsProp<T,radial>::doProperty(){
	static bool firstTime{true};
	vector<vector<int> > mCluster=this->Perco->getCluster();
	Matrix co=this->Mt.getCO();
	Matrix oc=this->Mt.getOC();
	vector<Dvect> xcm(mCluster.size(),Dvect{T{0.0}}); // Geometric center of the solute residues

	vector<Gyration<T> *> Rg=vector<Gyration<T>*>(mCluster.size());
	for(auto & ip: Rg)
		ip=new Gyration<T>();
	Gyration<T>::setTime(this->time_c);

	this->CalcGyro(this->mass,Rg);
	T r{-1};
	for(size_t o{0};o<Rg.size();o++){
		double Rh{5.0*Rg[o]->gRadg()/3.0};
		if(r < Rh) r=Rh;
		Dvect xcm0=Rg[o]->gXcm();
		xcm[o]=oc*xcm0;
	}
	if(firstTime){
		rcut=sqrt(r)*unit_nm;
		vector<T> tmp{co[XX][XX],co[YY][YY],co[ZZ][ZZ]};
		T rMax=*std::max_element(tmp.begin(),tmp.end());
		rcut=rcut > rMax?rMax:rcut;
		firstTime=false;
	}

	// Initialize histogram only once.
	static struct myAlloc{
		myAlloc(map<string,Properties::RhoHistogram *> & h, double r,vector<string> strs){
			for(auto str: strs)
			h[str]=new Properties::RhoHistogram(r);}
	} Once(histograms,rcut,*this->ResList0);

	// Obtain the mass of each molecule selected
	vector<int> * CMRes{nullptr};
	static struct myMasses{
		myMasses(vector<string> atres, vector<vector<int> > Sel, vector<double> mass, map<string,double> & Masses){
				for(size_t p{0};p<Sel.size();p++){
					double tmass{0};
					for(size_t q{0}; q< Sel[p].size();q++){
						size_t n=Sel[p][q];
						tmass+=mass[n];
					}
					Masses[atres[Sel[p][0]]]=tmass;
				}

		}
	}OnceAllMasses(this->atres,this->SelRes,this->mass,masses);

	// Find out if CMRes is null, if not null apply only to a single cluster, cluster 0.
	if(!this->MyRes.empty()){
		CMRes=new vector<int>{this->MyRes};

		//  Compute xcm[0] if there is a non null CMRes
		xcm[0]=Dvect{0};
		int NN{0};
		for(size_t o{0};o<CMRes->size();o++){

			auto p=CMRes->at(o);
			for(size_t q{0}; q< this->SelRes[p].size();q++){
				size_t n=this->SelRes[p][q];
				xcm[0]+=this->xa[n];
				NN++;
			}

		}
		xcm[0]/=static_cast<double>(NN);
	}
	for(size_t o{0};o<mCluster.size();o++){
		for(size_t p{0};p<this->SelRes.size();p++){
			Dvect xa{0};
			for(size_t q{0}; q< this->SelRes[p].size();q++){
				size_t n=this->SelRes[p][q];
				xa+=this->xa[n];
			}
			string Type=this->atres[this->SelRes[p][0]];
			xa/=(T) this->SelRes[p].size();
			Dvect xa_c=xa-xcm[o];
			for(size_t r{0};r<DIM;r++){
				xa_c[r]=xa_c[r]-rint(xa_c[r]);
			}
			Dvect xc=co*xa_c;
			T dist=sqrt(xc[XX]*xc[XX]+xc[YY]*xc[YY]+xc[ZZ]*xc[ZZ]);
			(*histograms[Type])(dist);
		}

	}
	for(auto it=histograms.begin();it != histograms.end();it++){
		(*it->second)++;
	}
	// Define the molecule mass for each molecule once and for all
	static struct myType{
		myType(map<string,Properties::RhoHistogram *> & hist, map<string,double>  & mass){
			for(auto it=hist.begin();it != hist.end();it++){
				it->second->setMass(mass[it->first]);

			}
		}
	} OnceMass(histograms,masses);
	return &histograms;
};

template <typename T>
void * AtomsProp<T,gyro>::doProperty(){
	vector<vector<int> > mCluster=this->Perco->getCluster();
	vector<Gyration<T> *> Rg=vector<Gyration<T>*>(mCluster.size());
	for(auto & ip: Rg)
		ip=new Gyration<T>();
	Gyration<T>::setTime(this->time_c);

	this->CalcGyro(this->mass,Rg);
	this->template Gyro<Enums::noJSON>();

	out=std::make_tuple(this->getRg_i(),this->gPerco());
	return &out;
}

template <typename T>
void * AtomsProp<T,gyroJ>::doProperty(){
	vector<vector<int> > mCluster=this->Perco->getCluster();
	vector<Gyration<T> *> Rg=vector<Gyration<T>*>(mCluster.size());
	for(auto & ip: Rg)
		ip=new Gyration<T>();
	Gyration<T>::setTime(this->time_c);

	this->CalcGyro(this->mass,Rg);
	this->template Gyro<Enums::JSON>();

	out=std::make_tuple(this->getRg_i(),this->gPerco());
	return &out;
}
template <typename T>
void * AtomsProp<T,pdbclust>::doProperty(){
	return this;
}
template <typename T>
void * AtomsProp<T,pdb>::doProperty(){
	return this;
}

template <typename T>
void * AtomsProp<T,pdbavg>::doProperty(){
	vector<vector<int> > mCluster=this->Perco->getCluster();
	vector<vector<int> > mAtoms=this->Perco->getAtoms();
	Matrix co=this->Mt.getCO();
	for(size_t o=0;o<mCluster.size();o++){
		for(size_t p=0;p<mCluster[o].size();p++){

		}
	}
	return this;
}

template <typename T>
void AtomsProp<T,pdb>::cPrint(ostream & fout){
	vector<Dvect> xc=vector<Dvect>(this->nr);
	for(int n=0;n<this->nr;n++)
		xc[n]=this->x[n];


	stringstream ss;
	ss<<this->time_c;
	fout << string("REMARK    Generated by trjProp ") <<endl;
	fout << string("REMARK    SIMULATION TIME = "+ ss.str()) <<endl;

	vector<T> Par=this->Mt.getParas();
	ss.str(string());
	ss<< setw(9) << setprecision(3) << fixed << Par[0]*10.0;
	ss<< setw(9) << setprecision(3) << fixed << Par[1]*10.0 ;
	ss<< setw(9) << setprecision(3) << fixed << Par[2]*10.0 ;
	ss<< setw(7) << setprecision(2) << fixed << Par[3] ;
	ss<< setw(7) << setprecision(2) << fixed << Par[4] ;
	ss<< setw(7) << setprecision(2) << fixed << Par[5] ;
	fout<< string("CRYST1"+ ss.str() + " P 1           1") <<endl;
	fout << "MODEL "<< fixed << setw(5) << right << this->cPrint_calls<<endl;
	for(size_t n{0};n<this->nr;n++){
		string tmp;
		ss.str(string());
		ss <<std::setw(5) << std::left << this->PDB[n].resn;
		string tmp2=this->PDB[n].first;
		tmp2.erase(26,1);
		tmp=tmp2.replace(21,5,ss.str());
		ss.str(string());
		ss << std::setw(5) << std::right << this->PDB[n].res+1;
		tmp.replace(21,5,ss.str());
		tmp.append(" 1.00");
		tmp.append(this->PDB[n].last);
		ss.str(string());
		ss << fixed << setw(8) << setprecision(3) << right << xc[n][XX]/unit_nm;;
		ss << fixed << setw(8) << setprecision(3) << right << xc[n][YY]/unit_nm;;
		ss << fixed << setw(8) << setprecision(3) << right << xc[n][ZZ]/unit_nm;;
		tmp.replace(30,24,ss.str());
		fout << tmp<<endl;

	}
	fout << string("ENDMDL") <<endl;
}

template class AtomsProp<float,radial>;
template class AtomsProp<double,radial>;
template class AtomsProp<float,gyro>;
template class AtomsProp<double,gyro>;
template class AtomsProp<float,gyroJ>;
template class AtomsProp<double,gyroJ>;
template class AtomsProp<float,pdb>;
template class AtomsProp<double,pdb>;
template class AtomsProp<float,pdbclust>;
template class AtomsProp<double,pdbclust>;
