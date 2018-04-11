/* CenterMass.h	-- Massimo Marchi
 *
 * $Header$
 * $Log$
 */
#ifndef _CENTERMASS_H
#define _CENTERMASS_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "MyUtilClass.h"
#include "Quaternions.hpp"
#include "WignerDMatrices.hpp"
#include <sstream>
using namespace std;
using namespace DVECT;
using namespace MATRIX;
using vectorvi=vector<vector<int> >;
using Quaternions::Quaternion;
using SphericalFunctions::WignerDMatrix;
using  vvector_i=vector<vector<int> >;

// diffk: do atomic diffusion
// diffq: do quaternion diffusion of detergent molecules
// diffMicelles: do rotational and translational diffusion of micelles
enum CenterMass_t {diff0, diffk, diffq, diffMicelles};
const double PI=3.14159265359;

template <typename T>
	class CenterMass{
protected:
	using Dvect=DDvect<T>;
	using Matrix=MMatrix<T>;
	using vectorvd=vector<vector<Dvect> >;
	using  vvvector_d=vector<vector<vector<Dvect> > >;
	using  vvector_d=vector<vector<Dvect> >;
  vector<Dvect> cm;
  vector<vector<Dvect> > cmm;
  Dvect Axis;
  Matrix CO,OC;
  static double Dt;
  static double time_c;
  vector<Quaternion> Q_t;
  vector<Quaternion> Qm;

// Will include the CO, OC, atomic, c.of.m. and quaternion trajectories
  vector<vector<Dvect> > cmt;
  vector<vector<Quaternion> > Q_tt;
  vector<vector<Quaternion> > Qmt;
  vector<Matrix> COt;
  vector<Matrix> OCt;
// end
  virtual void WriteIt(ostream &);
  virtual void ReadIt(ifstream &);
  virtual void SkipIt(ifstream &)=0;
  bool DoneOnce;
  Matrix MatrixfromQ(Quaternions::Quaternion );
  void RotateBack(vector<vector<Dvect> > & );
  virtual CenterMass * doClone() const=0;
public:
  CenterMass(): DoneOnce{false} {cm.clear();cmm.clear();};
  static void setTime(double);
  void setCOs(Matrix co, Matrix oc){ CO(co);OC(oc);}
  void setAxis(T x, T y, T z){Axis[XX]=x;Axis[YY]=y;Axis[ZZ]=z;}
  virtual void setup(size_t m){
	  cm=vector<Dvect>(m,{0.0,0.0,0.0});
	  cmm=vectorvd(m);
	  Q_t=vector<Quaternion>(m,Quaternion{1.0,0.0,0.0,0.0});
  }
  virtual void setup(size_t m, size_t n){
	  setup(m);
  }
  Dvect & getAxis(){return Axis;}
  static double getTime(){return time_c;};
  virtual void operator()(vector<Dvect> &, vector<vector<Dvect> >&,  vector<Quaternion> &,  vector<Quaternion> &);
  virtual void operator()(vector<Dvect> &, vector<Quaternion> &);
  virtual void operator()(vector<Dvect> &);
  virtual void Diffusion(ofstream & fout){};
  virtual void DiffQ(ofstream & fout,const int a, const int b, const int c){};


  size_t Size(){return cm.size();};
  Dvect & operator[](const int i) {return cm[i];};
  vector<Quaternion> & getQ(){return Q_t;}
  vector<Dvect> & getR(){return cm;};
  virtual ~CenterMass(){};
  virtual void RefCoord(rvec * y, vector<string> &y0,vector<vector<int> > & z, vvector_i & z0){};
  virtual void GofR(ofstream &)=0;
  virtual CenterMass * Clone() const {return doClone();}

  template <typename Q>
  friend ostream & operator<<(ostream &, CenterMass<Q> &);
  template <typename Q>
  friend ifstream & operator>>(ifstream &, CenterMass<Q> & );
  template <typename Q>
  friend ifstream & operator+=(ifstream &,CenterMass<Q> &);
};


#endif
