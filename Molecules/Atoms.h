/*
 * Atoms.h
 *
 *  Created on: May 20, 2011
 *      Author: marchi
 */

#ifndef ATOMS_H_
#define ATOMS_H_

#include "Metric.h"
#include "AtomIndex.h"
#include<iostream>
#include <cmath>
#include <vector>
#include "MyUtilClass.h"
#include "HeaderTrj.h"
#include "Topol.h"
#include "TopolPDB.h"
#include "Contacts.h"
#include "Gyration.h"
#include "GyrationJSON.h"

#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_seek.h"
#include "Percolation.h"
#include "PercolationJSON.h"
#include "myEnums.hpp"
#include <CenterMassBW3.h>
#include "FindCell.h"
#include "properties.hpp"
#include "RhoHistogram.h"
#include "PBCvect.h"

using namespace DVECT;
using namespace MATRIX;
using std::vector;
using namespace Enums;
using std::map;


template <typename T>
class brot{
	MMatrix<T> co;
public:

	const MMatrix<T> & operator()(const double a,const double b,const double c,
			const double alfa,const double beta,const double gamma){
	double degrad=M_PI/180.0;

	double ax=a;
    double alf=cos(degrad*alfa);
    double bet=cos(degrad*beta);
    double qt=sin(degrad*gamma);
    double gam=cos(degrad*gamma);
    double bx=b*gam;
    double by=b*qt;
    double cx=c*bet;
    double cy=c*(alf-bet*gam)/qt;
    double cz=sqrt(c*c-cx*cx-cy*cy);
    co[YY][XX]=0.0;
    co[ZZ][XX]=0.0;
    co[ZZ][YY]=0.0;
    co[XX][XX]=ax;
    co[XX][YY]=bx;
    co[XX][ZZ]=cx;
    co[YY][YY]=by;
    co[YY][ZZ]=cy;
    co[ZZ][ZZ]=cz;
    return co;
	}

};

template <typename T>
class Atoms {
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
protected:
	static int cPrint_calls;
	static int cPDBavg_calls;
	bool firsttime{true};
	const T HALF{0.5000-0.0001};
	static int calls;
	static vector<string> * ResList0;
	static vector<int> * ResIndx0;
	static string MIons;
	static string Mdetg;

	static map<string,double> MMass;
	static map<string,double> MMassNCH;
	vector<string> atres;
	vector<string> atname;

	int nr{0};
	int status{0};
	vector<Dvect> x;
	vector<Dvect> xa;
	Metric<T> Mt;
	vector<Dvect> x_avg;
	Matrix nCO;
	vector<Gyration<T> *> Rg_i;
	static int    step_c;
	static float  prec_c,time_c;
	vector<double> rd;
	vector<double> mass,massNCH;
	int Rg_count{0};
	CenterMass<T> * R_cmx{nullptr};
	vector<int> TypeNo;
	vector<vector<size_t>> AtsPerRes;


	virtual void ReadaStep(FstreamC * );
	virtual void ReadaStep(FstreamF * );
	virtual void ReadaStep(std::ifstream &);
	virtual void ReadaStepF(std::ifstream &);
	virtual void moveOffset(FstreamC *);
	virtual void moveOffset(FstreamF *);
	virtual void moveOffset(std::ifstream &);
	virtual void WriteaStep(FstreamC * );
	virtual void WriteaStep(FstreamF * );
	virtual void __ReconstructOneCluster(vector<bool> &);
	Percolation<T> * Perco{nullptr};
	void CalcGyro(vector<double> &,vector<Gyration<T> *> &);
	Dvect __FindCell(const vector<vector<int>> & ,const vector<vector<int>> & );
	virtual void cPrint(ostream &);
	void ndxPrint(ostream &);
	Topol_NS::TopolPDB PDB;
	vector<Topol_NS::PDBdata> PDBs;

	bool bndx{false};
	vector<vector<int> > SaxsSolute;
	vector<vector<int> > SelRes,myPadding;
	vector<int> MyRes;
public:
	Atoms(): R_cmx{new CenterMassBW3<T>}{};
	Atoms(const int);
	Atoms(const AtomIndex &);
	Atoms(const Atoms &);
	virtual ~Atoms();
	struct plane{
		Typedefs::real xc[4];
		int n;
		plane & operator=(Typedefs::rvec & xa){for(int o=0;o<DIM;o++) xc[o]=xa[o];return *this;};
		plane & operator=(Typedefs::real dd){xc[3]=dd;return *this;};
		plane & operator=(double dd){xc[3]=dd;return *this;};
		plane & operator=(int i){n=i;return *this;};
	} ax;

	virtual void doTest(){};
	virtual void * doProperty(){return nullptr;};
	virtual void printProperty(){};
	virtual void Reduce(Parallel::NewMPI *){};
	virtual myOptions getProperty(){return myOptions::noprop;}

	template <Enums::myWriteOptions OPT>
	void Gyro();
	vector<Gyration<T>*> & getRg_i(){return Rg_i;};
	const vector<int> & getTypeNo() const {return TypeNo;}
	vector<int> & getIndx(){return *ResIndx0;}
	void setrd(vector<double> & rdx) {rd=rdx;};
	void setrd(Topol_NS::Topol & y){rd=y.getrd();};

	void setDim(const int n);
	void setCoord(const Metric<T> &, const rvec *, const AtomIndex & );
	void setCoord(const Metric<T> &, const rvec *);
	void setMT(const Metric<T> &);
	void setMyRes(vector<int>);
	Atoms &  operator=(const Atoms &);
	Atoms & operator=(const T);
	Dvect & operator[](const int i){return x[i];};
	Atoms & operator()(const int);
	void Reconstruct(Contacts<T> *);
	void setNdx(bool b){bndx=b;}
	void Rot(const Matrix);
	int getNR()const{return nr;};
	const Metric<T> & getMt() const {return Mt;};
	Atoms COtoOC();
	Atoms OCtoCO();
	void doCOtoOC();
	void doOCtoCO();
	Atoms Shift(const rvec);
	Typedefs::real Dist(const int,const int);
	vector<Dvect> getX() const {return x;};
	vector<Dvect> getXA() const {return xa;};
	void setTopol(Topol_NS::Topol &);
	double getrd(int n){return rd[n];}
	CenterMass<T> & getCM() const {return *R_cmx;}
	virtual bool CenterAtoms();
	const vector<vector<int> > & getSaxsSolute() const {return SaxsSolute;}
	template <myAtoms MM>
	void InitSelection(vector<string> & y, Topol_NS::Topol & MyTop);
	void initLists(Topol_NS::TopolPDB &,vector<string> &);
	void CompCM();
	float getTime();
	template <Enums::myWriteOptions OPT>
	void SetupPercolate(Topol_NS::Topol &x);
	void SetupPercolate();

	int Percolate();
	int Percolate(double y){
		this->Perco->setRcut(y);
		return this->Percolate();
	}

	bool hasPerco() const {if(Perco) return true;return false;}
	Percolation<T> * gPerco(){return Perco;}

	void pdb(const vector<string> & c);
	vector<Dvect> getGC();
	vector<Comp> & getComp(){return Perco->getClustComp();}
	const vector<vector<int> > & getCluster() const;
	const vector<vector<int> > & getAtoms() const;
	void PDBavg();
	void printPDBavg(ostream & );

	friend Fstream & operator+=(Fstream & fin, Atoms & y){
		if(FstreamC * finC=dynamic_cast<FstreamC *> (&fin))
			y.moveOffset(finC);
		else if(FstreamF * finF=dynamic_cast<FstreamF *> (&fin))
			y.moveOffset(finF);

		return fin;
	}

	friend std::ifstream & operator+=(std::ifstream & fin,Atoms & y){
		y.moveOffset(fin);
		return fin;
	}

	friend Fstream & operator>>(Fstream & fin, Atoms & y){
		if(FstreamC * finC=dynamic_cast<FstreamC *> (&fin))
			y.ReadaStep(finC);
		else if(FstreamF * finF=dynamic_cast<FstreamF *> (&fin))
			y.ReadaStep(finF);
		return fin;
	}
	friend Fstream & operator<<(Fstream & fout, Atoms & y){
		if(FstreamC * foutC=dynamic_cast<FstreamC *> (&fout))
			y.WriteaStep(foutC);
		return fout;
	}
	friend ostream & operator<<(ostream & fout , Atoms & y ){
		y.cPrint(fout);
		return fout;
	}
	friend std::ifstream & operator>>(std::ifstream & fin, Atoms & y){
		y.ReadaStep(fin);
		return fin;
	}

};

#endif /* ATOMS_H_ */

