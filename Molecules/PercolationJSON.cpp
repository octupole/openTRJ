/*
 * PercolationJSON.cpp
 *
 *  Created on: Dec 21, 2017
 *      Author: marchi
 */

#include "PercolationJSON.h"

template <typename T>
PercolationJSON<T>::~PercolationJSON() {
	// TODO Auto-generated destructor stub
}

template <typename T>
void PercolationJSON<T>::__Writeit(ostream & fout){
	map<string,int> mapRes;
	map<string,int> mapResAtm;
	std::hash<string> str_hash;

	for(size_t o=0;o<this->Clusters.size();o++){
		json myClust;
		mapRes.clear();mapResAtm.clear();
		map<string,vector<int>> Tag;
		for(size_t p=0;p<this->Clusters[o].size();p++){
			int n=this->Clusters[o][p];
			mapRes[this->pResn[n]]++;
			Tag[this->pResn[n]].push_back(n);
			for(size_t q{0};q<this->Atoms[n].size();q++){
				mapResAtm[this->pResn[n]]++;
			}
		}
		json tmp;
		for(auto it{mapRes.begin()};it!= mapRes.end();it++){
			tmp[it->first]={it->second,mapResAtm[it->first]};
		}
		for(auto it{Tag.begin()};it!= Tag.end();it++){
			std::stringstream ss;
			vector<int> v=it->second;
			std::sort(v.begin(),v.end());
			std::copy(v.begin(),v.end(),std::ostream_iterator<int>( ss," "));
			tmp["hsh"][it->first]=str_hash(ss.str());
		}
		myJson["cluster"].push_back(tmp);
	}
	vector<T> co,oc;
	for(size_t n{0};n< DIM;n++)
		for(size_t m{0};m< DIM;m++){
			co.push_back(this->CO[n][m]);
			oc.push_back(this->OC[n][m]);
		}
	myJson["CO"]=co;
	myJson["OC"]=oc;
	fout <<"{";
	for(auto it=myJson.begin();it != myJson.end();++it){
		fout << "\""<<it.key() <<"\": ";
		fout<<myJson[it.key()]<<",";
	}
	myJson.clear();
}

template class PercolationJSON<float>;
template class PercolationJSON<double>;
