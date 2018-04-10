/*
 * CenterMassMicelles.h
 *
 *  Created on: Jun 3, 2015
 *      Author: marchi
 */

#ifndef SRC_CENTERMASSMICELLES_H_
#define SRC_CENTERMASSMICELLES_H_

#include "CenterMassBW3.h"
#include "histograms.hpp"

struct Histogram1Db: public Histogram1D{
	Histogram1Db(): Histogram1D(){};
	Histogram1Db(double dx0,double cut0);
	virtual hist1D & operator[](size_t );
	void Normalize(double);
	size_t size();
	friend ostream & operator<<(ostream &, Histogram1Db &);
	virtual ~Histogram1Db(){};

};

template <typename T>
class CenterMassMicelles: public CenterMassBW3<T>{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	using vectorvd=vector<vector<Dvect> >;
	using  vvvector_d=vector<vector<vector<Dvect> > >;
	using  vvector_d=vector<vector<Dvect> >;
	using CenterMassBW3<T>::DiffQQ;
	void CorrectCoords(vvector_d &);
	Dvect Rint(Dvect v){
		for(auto o=0;o<DIM;o++)
			v[o]=rint(v[o]);
		return v;
	}
public:
	CenterMassMicelles();
	virtual void DiffQ(ofstream &,const int , const int , const int);
	virtual void Diffusion(ofstream & fout);
	virtual void GofR(ofstream &);
	virtual ~CenterMassMicelles();
};

#endif /* SRC_CENTERMASSMICELLES_H_ */
