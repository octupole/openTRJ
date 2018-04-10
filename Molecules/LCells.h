/*
 * LCells.h
 *
 *  Created on: Dec 18, 2014
 *      Author: marchi
 */

#ifndef LINKEDCELLS_SRC_LCELLS_H_
#define LINKEDCELLS_SRC_LCELLS_H_

#include <vector>
#include <algorithm>
#include "Ftypedefs.h"
#include "MyUtilClass.h"



using namespace std;
using namespace MATRIX;
using namespace DVECT;


template <typename T>
class LCells {
	const T HALF{0.5000-0.0001};
	const T small{1.0e-6};
	using Matrix=MMatrix<T>;
	using Dvect=DDvect<T>;
	typedef vector<int> vectint;
	int nr{0};
	T Rcut{1.3};
	const T Rmax{4.5};

	vector<int> nc={-1,-1,-1};
	vector<vectint>  indx;
	vector<vectint> Cellp;
	vector<vector<vector<vectint> > > Chainp;
	Matrix co,oc;
	vector<Dvect> x;
	vector<vector<int> > nnl;
	T Dist_ijk(int, int,int);
	void Init(Matrix & co0, const vector<Dvect> & y);
public:
	LCells();
	LCells(Matrix & co0, const vector<Dvect> & y, T rcut): Rcut{rcut}{
		this->Init(co0,y);
		this->Index();
	}
	void operator()(Matrix & co0, const vector<Dvect> y){
		this->Init(co0,y);
	}
	vector<vector<int> > & getNeigh(){
		try{
			if(!nnl.size()) throw " Neighbour list not yet computed !";
		}
		catch(const char * s){
			cout << s << endl;
			exit(1);
		}
		return nnl;
	}
	void Index();
	bool test();
	vector<vector<int> > &List(bool=true);
	virtual ~LCells();
};

#endif /* LINKEDCELLS_SRC_LCELLS_H_ */
