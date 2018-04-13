/*
 * FindCell.cpp
 *
 *  Created on: May 18, 2016
 *      Author: marchi
 */

#include "FindCell.h"

template <typename T>
DDvect<T> FindCell<T>::Run(const vector<vector<int> > & mCluster, const vector<vector<int> > & mAtoms,const vector<Dvect> xa, int na){
	vector<Dvect> x(na);
	for(auto o=0;o<na;o++){
		for(auto n=0;n<DIM;n++){
			x[o][n]=xa[o][n];
		}
	}
	vector<T> xcm(mAtoms.size(),0.0);
	Dvect myReturn{T{0.0}};
	T dx=T{0.20};
	for(auto n=0;n<DIM;n++){
		size_t Nn=CO[n][n]/dx;
		size_t M{0};
		T ddx=1.0/Nn;
		T TT{T{0}};
		vector<int> myPBC(Nn);
		while(M < Nn){
			size_t bCount{0};
			for(auto o=0;o<na;o++){
				x[o][n]+=ddx;
			}
			for(size_t o=0;o<mAtoms.size();o++){
				// Compute molecule center of mass xcm
				xcm[o]=0.0;
				for(size_t p=0;p<mAtoms[o].size();p++){
					double xx=x[mAtoms[o][p]][n];
					xcm[o]+=xx;
				}
				xcm[o]/=static_cast<double>(mAtoms[o].size());
				xcm[o]-=rint(xcm[o]-0.5);
			}

			for(size_t o=0;o<mCluster.size();o++){
				int nn=mCluster[o][0];
				T Xref=xcm[nn];
				for(size_t p=1;p<mCluster[o].size();p++){
					int nn=mCluster[o][p];
					T PBC=rint(Xref-xcm[nn]);
					if(PBC) bCount++;
				}
			}
			myPBC[M]=bCount;
			M++;
			TT+=ddx;
		}
		auto result = std::minmax_element(myPBC.begin(),myPBC.end());
		auto it1=std::find_if(myPBC.begin(),myPBC.end(),[result](int q){return q==*result.first;});
		auto it2=std::find_if(it1,myPBC.end(),[result](int q){return q!=*result.first;});
		auto j=std::distance(myPBC.begin(),it1)+std::distance(it1,it2)/2;
		myReturn[n]=j*ddx;
	}
	return 	myReturn;
}
template <typename T>
FindCell<T>::~FindCell() {
	// TODO Auto-generated destructor stub
}
template class FindCell<float>;
template class FindCell<double>;
