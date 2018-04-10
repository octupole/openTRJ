/*
 * Contacts.h
 *
 *  Created on: May 28, 2015
 *      Author: marchi
 */

#ifndef SRC_CONTACTS_H_
#define SRC_CONTACTS_H_
#include "LCells.h"
#include <algorithm>
#include <limits>
const size_t SIZE_T=10000000;
#define RCUT 1.5
#define RCUT_IN 1.3
template <typename T>
class Contacts {
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
protected:
	vector<Dvect> v;
	vector<bool> b;
	vector<size_t> List;
	size_t Ind=0;
	vector<vector<int> > nnl;

	Matrix CO,OC;
	T Rcut{RCUT}, Rcut_in{RCUT_IN};
	virtual void Init(vector<Dvect> & y,Matrix & co, Matrix & oc){
		CO=co;OC=oc;v.clear();v=y;b.clear();b=vector<bool>(v.size(),true);
		List.clear();
		Ind=0;
	}
	void rCluster(size_t);
	void CompNei();
public:
	Contacts(vector<Dvect> & y,Matrix & co, Matrix & oc, T R, T Ri, size_t mystart): Ind(mystart), Rcut(R*0.1), Rcut_in(Ri*0.1){
		Init(y,co,oc);
	}
	Contacts(vector<Dvect> & y,Matrix & co, Matrix & oc, T R, T Ri): Rcut(R*0.1), Rcut_in(Ri*0.1) {
		Init(y,co,oc);
	}
	Contacts(vector<Dvect> & y,Matrix & co, Matrix & oc, T R): Rcut(R*0.1) {
		Init(y,co,oc);
	}
	Contacts(vector<Dvect> & y,Matrix & co, Matrix & oc){
		Init(y,co,oc);
	}
	Contacts(T R,T R0):  Rcut(R*0.1),Rcut_in(R0*0.1){}
	Contacts(T R0):  Rcut_in(R0*0.1){}
	Contacts(){};
	void operator()(vector<Dvect> & y,Matrix & co, Matrix & oc){
		Init(y,co,oc);
	}
	vector<vector<int>> & NNL(){return nnl;}
	void setR(T R, T Ri){Rcut=R*0.1;Rcut_in=Ri*0.1;}
	virtual void Neighbors();
	size_t next();
	Dvect & operator[](size_t n){return v[n];}
	virtual ~Contacts();
};
#endif /* SRC_CONTACTS_H_ */
