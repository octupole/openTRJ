/*
 * PBCvect.h
 *
 *  Created on: Mar 10, 2019
 *      Author: marchi
 */

#ifndef MOLECULES_PBCVECT_H_
#define MOLECULES_PBCVECT_H_
template <typename T>
class PBCvect{
	using Dvect=DDvect<T>;
	static vector<Dvect> * nxyz;
	PBCvect(){};
	static void __compVect(){
		const int M{3};
		nxyz=new vector<Dvect>;
		for(int o=0;o<M;o++){
			double nx=static_cast<double>(o-1);
			for(int p=0;p<M;p++){
				double ny=static_cast<double>(p-1);
				for(int q=0;q<M;q++){
					double nz=static_cast<double>(q-1);
					nxyz->push_back(Dvect(nx,ny,nz));
				}
			}
		}

	};
public:
	static vector<Dvect> & getVec(){__compVect();return *nxyz;}
};

template <typename T>
vector<DDvect<T>> * PBCvect<T>::nxyz=nullptr;





#endif /* MOLECULES_PBCVECT_H_ */
