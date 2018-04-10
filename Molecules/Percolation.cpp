/*
 * Percolation.cpp
 *
 *  Created on: Sep 12, 2013
 *      Author: marchi
 */

#include "Percolation.h"


template <typename T>
vector<string> Percolation<T>::Select=vector<string>{"S","OS1","OS2","OS3","NA"};

template <typename T>
double Percolation<T>::PercoCutoff{0.0};

const double OFFSET=1.5;

struct op_mysort{
	bool operator()(const vector<int>  & x,const vector<int> & y){
		return x.size() > y.size();
	}
};

template <typename T>
Percolation<T>::Percolation(const listcon & y,const  vector<string> & atnme, const vector<double>  & rd){
	string tmp;
	int m=0;
	do{tmp+=" "+Select[m++]+" ";}while(m < Select.size());

	nr=0;
	for(size_t o=0; o < y.size() ;o++){
		vector<int> atms;
		vector<string> atname;
		for(size_t p=0; p <y[o].size() ; p++){
			int n=y[o][p];
			if(tmp.find(atnme[n]) != string::npos){
				atms.push_back(n);
				atname.push_back(atnme[n]);
			}
		}
		if(atms.empty()) continue;
		if(atms.size() == 1) Name.push_back("NA");
		else Name.push_back("AOT");
		Atoms.push_back(atms);
		pAtn.push_back(atname);
		nr++;
	}
	pRd=rd;
	myPercoCutoff=[this](int i,int j){if(PercoCutoff){return PercoCutoff;};return (pRd[i]+pRd[j])*OFFSET;};
}
template <typename T>
Percolation<T>::Percolation(const vector<vector<int> >& y,const  vector<double>  & rd,const   vector<string> & resn){
	string tmp;
	int m=0;
	nr=rd.size();
	pRd=rd;
	Atoms=y;

	pResn=vector<string>(y.size());
	for(size_t o=0;o<Atoms.size();o++){
		pResn[o]=resn[Atoms[o][0]];
	}
	myPercoCutoff=[this](int i,int j){if(PercoCutoff){return PercoCutoff;};return (pRd[i]+pRd[j])*OFFSET;};
}
template <typename T>
Percolation<T>::Percolation(const vector<vector<int> >& y,const  vector<double>  & rd,const   vector<string> & resn, const vector<string> & atms) {
	string tmp;
	int m=0;
	nr=rd.size();
	pRd=rd;
	Atoms=y;
	pResn=vector<string>(y.size());
	for(size_t o=0;o<Atoms.size();o++){
		pResn[o]=resn[Atoms[o][0]];
	}
	pAtn=vector<vector<string> >(y.size());

		int nato{0};
	for(size_t o=0;o<Atoms.size();o++){
		pAtn[o]=vector<string>(Atoms[o].size());
		for(size_t p=0;p<Atoms[o].size();p++){
			int q=Atoms[o][p];
			pAtn[o][p]=atms[q];
		}
	}
	myPercoCutoff=[this](int i,int j){if(PercoCutoff){return PercoCutoff;};return (pRd[i]+pRd[j])*OFFSET;};
}

template <typename T>
void Percolation<T>::rCluster(size_t m){
	bAtoms[m]=m;
	for(size_t n=0; n < Contacts[m].size(); n++){
		size_t na=Contacts[m][n];
		if(bAtoms[na] < 0) rCluster(na);
	}
	return;
}
template <typename T>
int  Percolation<T>::gCluster(){
	bAtoms.clear();
	Clusters.clear();
	bAtoms=vector<int>(Atoms.size(),-1);
	vector<int> bbAtoms=vector<int>(Atoms.size());
	size_t n=0;
	int ma=0;
	do{
		bbAtoms=bAtoms;
		vector<int> cluster;
		if(bAtoms[n] < 0) rCluster(n);

		vector<int> tmp1;
		tmp1.assign(bAtoms.begin(),bAtoms.end());
		vector<int> tmp2;
		tmp2.assign(bbAtoms.begin(),bbAtoms.end());
		std::sort(tmp1.begin(),tmp1.end());
		std::sort(tmp2.begin(),tmp2.end());
		std::set_difference(tmp1.begin(),tmp1.end(),tmp2.begin(),tmp2.end(),std::back_inserter(cluster));
		Clusters.push_back(cluster);
		do{
			n++;
		} while(bAtoms[n] > -1 && n < bAtoms.size() );
		ma++;
	} while(n < bAtoms.size());
	std::sort(Clusters.begin(),Clusters.end(),op_mysort());
	ClustComp.clear();
	ClustComp=vector<Comp>(Clusters.size());
	for(size_t o=0;o<Clusters.size();o++){
		map<string,vector<int> > maps;
		for(size_t p=0;p<Clusters[o].size();p++){
			maps[pResn[Clusters[o][p]]].push_back(Clusters[o][p]);
		}
		map<string,vector<int> >::iterator ip=maps.begin();
		for(;ip!=maps.end();ip++){
			int n=ip->second.size();
			ClustComp[o].Res.push_back(ip->first);
			ClustComp[o].No.push_back(n);
		}

	}
	return Clusters.size();
}

template <typename T>
void Percolation<T>::doContacts(vector<Dvect> & v){
	Contacts.clear();

	Contacts=vector<vector<int> >(nr);

	for(size_t n=0; n < Atoms.size() ; n++){
		vector<int> tmp;
		for(size_t o=0; o < Atoms[n].size(); o++){
			int i=Atoms[n][o];
			bool isNumeric=pAtn[n][o].substr(0,1) == "1" ||pAtn[n][o].substr(0,1) == "2" ||pAtn[n][o].substr(0,1) == "3";
			if(pAtn[n][o] != "C" && exclusion.find(pAtn[n][o]) != string::npos) continue;
			if(pAtn[n][o].substr(0,1) =="H" || (isNumeric && pAtn[n][o].substr(1,1) == "H")) continue;

			Dvect xm=v[i];
			for(size_t m=n+1; m < Atoms.size(); m++){
				for(size_t p=0; p < Atoms[m].size();p++){
					int j=Atoms[m][p];
					bool found2=exclusion.find(pAtn[m][p]) != string::npos;
					if(found2) continue;
					if(pAtn[m][p].substr(0,1) =="H") continue;

					Dvect xn=v[j];
					double dist=xm.Dist(xn);
					if(dist <= myPercoCutoff(i,j)){
						tmp.push_back(m);
					}
				}
			}
		}
		if(!tmp.empty()){
			std::sort(tmp.begin(),tmp.end());
			vector<int>::iterator it=std::unique(tmp.begin(),tmp.end());
			tmp.resize(std::distance(tmp.begin(),it));
			Contacts[n]=tmp;
		}
	}
	for(size_t n=0; n < Contacts.size(); n++){
		for(size_t m=0; m < Contacts[n].size(); m++) {
			if(Contacts[n][m] > n) Contacts[Contacts[n][m]].push_back(n);
		}

	}
}
template <typename T>
void Percolation<T>::doContacts(vector<Dvect> & v, Matrix & co, Matrix & oc){
	Contacts.clear();
	Contacts=vector<vector<int> >(nr);
	vector<Dvect> X;
	vector<vector<int> > N;
	for(size_t n=0; n < Atoms.size() ; n++){
		vector<int> tmp;
		for(size_t o=0; o < Atoms[n].size(); o++){
			int i=Atoms[n][o];
			bool isNumeric=pAtn[n][o].substr(0,1) == "1" ||pAtn[n][o].substr(0,1) == "2" ||pAtn[n][o].substr(0,1) == "3";
			if(pAtn[n][o] != "C" && exclusion.find(pAtn[n][o]) != string::npos) continue;
			if(pAtn[n][o].substr(0,1) =="H" || (isNumeric && pAtn[n][o].substr(1,1) == "H")) continue;
			Dvect xm=oc*v[i];
			X.push_back(xm);
			vector<int> tm=vector<int>(2);
			tm[0]=n;
			tm[1]=o;
			N.push_back(tm);
		}
	}
	LCells<T> Cells(co,X,Rcut);
	Cells.Index();
	vector<vector<int> > & nnl=Cells.List(false);
	vector<vector<int> > tmp=vector<vector<int> >(Atoms.size());

	for(size_t nn=0;nn<nnl.size();nn++){
		int n=N[nn][0];
		int o=N[nn][1];
		int i=Atoms[n][o];
		Dvect xm=co*X[nn];

		for(size_t mm0=0;mm0<nnl[nn].size();mm0++){
			int mm=nnl[nn][mm0];
			int m=N[mm][0];
			int p=N[mm][1];
			if(m == n) continue;
			int j=Atoms[m][p];
			Dvect xn=co*X[mm];
			double dist=xm.Dist(xn,co,oc);
			if(dist <= myPercoCutoff(i,j)){
				tmp[n].push_back(m);
			}
		}
	}
	for(size_t n=0; n < Atoms.size() ; n++){
		if(!tmp[n].empty()){
			std::sort(tmp[n].begin(),tmp[n].end());
			vector<int>::iterator it=std::unique(tmp[n].begin(),tmp[n].end());
			tmp[n].resize(std::distance(tmp[n].begin(),it));
			Contacts[n]=tmp[n];
		}
	}
	CO=co;
	OC=oc;
}
template <typename T>
void Percolation<T>::__Writeit(ostream & fout){
	map<string,int> mapRes;
	map<string,int> mapResAtm;

	for(size_t o=0;o<Clusters.size();o++){
		mapRes.clear();mapResAtm.clear();
		for(size_t p=0;p<this->Clusters[o].size();p++){
			int n=this->Clusters[o][p];
			mapRes[this->pResn[n]]++;
			for(size_t q{0};q<this->Atoms[n].size();q++){
				mapResAtm[this->pResn[n]]++;
			}
		}
		int ntot{0};
		for(auto it{mapRes.begin()};it!= mapRes.end();it++){
			ntot+=it->second;
		}
		fout << "  Cluster = "<< std::fixed << std::setw(3) << o << " Size " << ntot << " : ";
		for(auto it{mapRes.begin()};it!= mapRes.end();it++){
			fout << it->first << "[" << it->second << "] ";
		}
		fout << endl;
	}
}



template <typename T>
Percolation<T>::~Percolation() {
	// TODO Auto-generated destructor stub
}

template class Percolation<float>;
template class Percolation<double>;
