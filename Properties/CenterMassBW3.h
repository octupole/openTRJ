/*
 * CenterMassBW3.h
 *
 *  Created on: May 12, 2015
 *      Author: marchi
 */

#ifndef SRC_CENTERMASSBW3_H_
#define SRC_CENTERMASSBW3_H_

#include "CenterMass.h"
#include <linalg.h>
#include "fasttransforms.h"
#include "CMPick.h"
#include <map>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include "Timer.h"
using vvector_q=vector<vector<Quaternion> >;

struct mysort{
	bool operator()(const string & i, const string & j){
		if(i[0] != j[0])
			return i>j;
		else{
			if(i.size() != j.size())
				return i.size() <j.size() ;
			else
				return i<j;
		}
	}
};
template <typename T>
struct dual_d{
	T x,y;
};
template <typename T>
class CenterMassBW3: public CenterMass<T>{
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	using vectorvd=vector<vector<Dvect> >;
	using  vvvector_d=vector<vector<vector<Dvect> > >;
	using  vvector_d=vector<vector<Dvect> >;
	using CenterMass<T>::time_c;
protected:
	  virtual void WriteIt(ostream &);
	  virtual void ReadIt(ifstream &);
	  virtual void SkipIt(ifstream &);

	  void RotateTranslate(vvector_d &, vvector_d & );
	  static string HeavyAtoms1;
	  static string HeavyAtoms2;
	  vector<string> Labels;
	  vvector_d xc; // Rotation corrected atomic coordinates
	  vvector_d xd; // Atomic rotations
	  vvector_d X;
	  vvector_d Xd; // Atomic rotations read in

	  map<string,vector<dual_d<T>>, mysort> Diffs;
	  map<string,long int, mysort> nDiffs;
	  static bool BeenThere;
	  alglib::real_1d_array Diff1D(vector<Dvect> & , size_t );
	  virtual vector<complex<T>> DiffQQ(vvector_q &,const int , const int , const int);
	  const T R00=1.3;
	  bool AlreadySkipped=false;
	  streampos SkipPos=0;
	  virtual CenterMassBW3 * doClone() const {return new CenterMassBW3(*this);};

public:
	CenterMassBW3();
	virtual ~CenterMassBW3();
	virtual void setup(size_t m, size_t n){
		  CenterMass<T>::setup(m);
		  this->Qm=vector<Quaternion>(n,Quaternion{1.0,0.0,0.0,0.0});
	  }
	virtual CenterMassBW3 * Clone() const {return doClone();}

	virtual void operator()(vector<Dvect> &, vector<vector<Dvect> >&,  vector<Quaternion> &, vector<Quaternion> &);
	virtual void RefCoord(rvec *, vector<string> &,vector<vector<int> > &, vvector_i & );
	virtual void Diffusion(ofstream & fout);
	virtual void DiffQ(ofstream & fout,const int a, const int b, const int c){};
	virtual void GofR(ofstream & fout){};
};


#endif /* SRC_CENTERMASSBW3_H_ */
