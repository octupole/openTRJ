/*
 * AtomsProp.cpp
 *
 *  Created on: Mar 1, 2019
 *      Author: marchi
 */

#include "AtomsProp.h"
template <typename T>
string AtomsProp<T,radial>::headerXVG=R"(
# Grace project file
#
@version 50122
@page size 792, 612
@page scroll 5%
@page inout 5%
@link page off
@map font 13 to "ZapfDingbats", "ZapfDingbats"
@map font 4 to "Helvetica", "Helvetica"
@map font 6 to "Helvetica-Bold", "Helvetica-Bold"
@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"
@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"
@map font 17 to "Nimbus-Sans-L-Condensed-Regular", "Nimbus-Sans-L-Condensed-Regular"
@map font 18 to "Nimbus-Sans-L-Condensed-Bold", "Nimbus-Sans-L-Condensed-Bold"
@map font 19 to "Nimbus-Sans-L-Condensed-Regular-Italic", "Nimbus-Sans-L-Condensed-Regular-Italic"
@map font 20 to "Nimbus-Sans-L-Condensed-Bold-Italic", "Nimbus-Sans-L-Condensed-Bold-Italic"
@map font 0 to "Times-Roman", "Times-Roman"
@map font 2 to "Times-Bold", "Times-Bold"
@map font 1 to "Times-Italic", "Times-Italic"
@map font 3 to "Times-BoldItalic", "Times-BoldItalic"
@map font 8 to "Courier", "Courier"
@map font 10 to "Courier-Bold", "Courier-Bold"
@map font 9 to "Courier-Oblique", "Courier-Oblique"
@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"
@map font 29 to "URW-Palladio-L-Roman", "URW-Palladio-L-Roman"
@map font 30 to "URW-Palladio-L-Bold", "URW-Palladio-L-Bold"
@map font 31 to "URW-Palladio-L-Italic", "URW-Palladio-L-Italic"
@map font 32 to "URW-Palladio-L-Bold-Italic", "URW-Palladio-L-Bold-Italic"
@map font 12 to "Symbol", "Symbol"
@map font 34 to "URW-Chancery-L-Medium-Italic", "URW-Chancery-L-Medium-Italic"
@map color 0 to (255, 255, 255), "white"
@map color 1 to (0, 0, 0), "black"
@map color 2 to (255, 0, 0), "red"
@map color 3 to (0, 255, 0), "green"
@map color 4 to (0, 0, 255), "blue"
@map color 5 to (255, 255, 0), "yellow"
@map color 6 to (188, 143, 143), "brown"
@map color 7 to (220, 220, 220), "grey"
@map color 8 to (148, 0, 211), "violet"
@map color 9 to (0, 255, 255), "cyan"
@map color 10 to (255, 0, 255), "magenta"
@map color 11 to (255, 165, 0), "orange"
@map color 12 to (114, 33, 188), "indigo"
@map color 13 to (103, 7, 72), "maroon"
@map color 14 to (64, 224, 208), "turquoise"
@map color 15 to (0, 139, 0), "green4"
@reference date 0
@date wrap off
@date wrap year 1950
@default linewidth 1.0
@default linestyle 1
@default color 1
@default pattern 1
@default font 0
@default char size 1.000000
@default symbol size 1.000000
@default sformat "%.8g"
@background color 0
@page background fill on
@timestamp off
@timestamp 0.03, 0.03
@timestamp color 1
@timestamp rot 0
@timestamp font 0
@timestamp char size 1.000000
@timestamp def "Mon Oct  8 15:33:20 2012"
)";


template <typename T>
void AtomsProp<T,radial>::Reduce(Parallel::NewMPI * y){
	for(auto it{histograms.begin()};it != histograms.end();it++){
		it->second->Reduce(y);
	}
}


template <typename T>
void * AtomsProp<T,radial>::doProperty(){
	static bool firstTime{true};
	vector<vector<int> > mCluster=this->Perco->getCluster();
	Matrix co=this->Mt.getCO();
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
		xcm[o]=Rg[o]->gXcm();
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
	for(auto it{histograms.begin()};it != histograms.end();it++){
		(*it->second)++;
	}
	// Define the molecule mass for each molecule once and for all
	static struct myType{
		myType(map<string,Properties::RhoHistogram *> & hist, map<string,double>  & mass){
			for(auto it{hist.begin()};it != hist.end();it++){
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
void * AtomsProp<T,pdb>::doProperty(){
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

template class AtomsProp<float,radial>;
template class AtomsProp<double,radial>;
template class AtomsProp<float,gyro>;
template class AtomsProp<double,gyro>;
template class AtomsProp<float,gyroJ>;
template class AtomsProp<double,gyroJ>;
