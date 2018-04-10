/*
 * Contacts.cpp
 *
 *  Created on: May 28, 2015
 *      Author: marchi
 */

#include "Contacts.h"

template <typename T>
void Contacts<T>::rCluster(size_t m){
	b[m]=false;
	for(size_t n=0; n < nnl[m].size(); n++){
		size_t na=nnl[m][n];
		if(b[na]) {
			List.push_back(na);
			rCluster(na);
		}
	}
	return;
}
template <typename T>
void Contacts<T>::CompNei(){
	T cut2=Rcut_in*Rcut_in;
	nnl=vector<vector<int>>(v.size());
	for(size_t o=0;o<v.size();o++){
		Dvect x_o=v[o];
		for(size_t p=0;p<v.size();p++){
			if(p == o) continue;
			Dvect x_p=v[p];
			Dvect xs{x_p-x_o};
			for(int q=0;q<DIM;q++){
				xs[q]-=rint(xs[q]);
			}
			Dvect xc=CO*xs;
			if(xc[XX]*xc[XX]+xc[YY]*xc[YY]+xc[ZZ]*xc[ZZ] < cut2){
				nnl[o].push_back(p);
			}
		}
	}
}
template <typename T>
void Contacts<T>::Neighbors(){
	List.push_back(Ind);
	LCells<T> Cells(CO,v,Rcut);

	if(v.size() > 1){
		if(Cells.test()){
			Cells.Index();
			nnl=Cells.List(false);
		} else{
			CompNei();
		}

		rCluster(Ind);
	}
}
template <typename T>
size_t Contacts<T>::next(){
	if(Ind < List.size()){
		size_t i_out=List[Ind];
		Ind++;
		return i_out;
	}
	else
		return SIZE_T;
}

template <typename T>
Contacts<T>::~Contacts() {
	// TODO Auto-generated destructor stub
}
template class Contacts<float>;
template class Contacts<double>;
