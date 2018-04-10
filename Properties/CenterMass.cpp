
#include "CenterMass.h"

template <typename T>
double CenterMass<T>::time_c=-1.0;
template <typename T>
double CenterMass<T>::Dt=0.0;

template <typename T>
void CenterMass<T>::WriteIt(ostream & fout) {
	  size_t beg0=0, end0=cm.size();
	 // if( rcm.Size() != ) return fout;
	  fout << fixed << setw(8) << setprecision(2) << CenterMass::getTime() << " : ";
	  for(size_t o=beg0;o<end0;o++){
	    fout << fixed << setw(9) << setprecision(5) << cm[o][XX] << " ";
	    fout << fixed << setw(9) << setprecision(5) << cm[o][YY] << " ";
	    fout << fixed << setw(9) << setprecision(5) << cm[o][ZZ] << " : ";
	  }
	  fout << Axis[XX] << " " << Axis[YY] << " " << Axis[ZZ] ;
	  fout << endl;

}
template <typename T>
void CenterMass<T>::ReadIt(ifstream & fin){

}
template <typename T>
ifstream & operator+=(ifstream & fin,CenterMass<T> & rcm){
	rcm.SkipIt(fin);
	return fin;
}
template <typename T>
ifstream & operator>>(ifstream & fin, CenterMass<T> & rcm){
	rcm.ReadIt(fin);
	return fin;
}
template <typename T>
ostream & operator<<(ostream & fout , CenterMass<T> & rcm){
	rcm.WriteIt(fout);
	return fout;
}
template <typename T>
MMatrix<T> CenterMass<T>::MatrixfromQ(Quaternions::Quaternion Q){
	Matrix mat;
	double W=Q[0],X=Q[1],Y=Q[2],Z=Q[3];
	double xx,xy,xz,xw,yy,yz,yw,zz,zw;
	xx=X*X;
	xy=X*Y;
	xz=X*Z;
	xw=X*W;
	yy=Y*Y;
	yz=Y*Z;
	yw=Y*W;
	zz=Z*Z;
	zw=Z*W;

	mat[XX][XX]=1-2*(yy+zz);
	mat[XX][YY]=2*(xy-zw);
	mat[XX][ZZ]=2*(xz+yw);
	mat[YY][XX]=2*(xy+zw);
	mat[YY][YY]=1-2*(xx+zz);
	mat[YY][ZZ]=2*(yz-xw);
	mat[ZZ][XX]=2*(xz-yw);
	mat[ZZ][YY]=2*(yz+xw);
	mat[ZZ][ZZ]=1-2*(xx+yy);
	return mat;
}
template <typename T>
void CenterMass<T>::RotateBack(vector<vector<Dvect> > & X) {
	Matrix Mat;
	for(auto p=0;p<X.size();p++){
		Quaternion qq=Q_t[p];
		Matrix Rot=MatrixfromQ(qq.inverse());
		for(auto q=0;q<X[p].size();q++){
			Dvect tmp=Rot*X[p][q];
			X[p][q]=tmp;
		}
	}
}

template <typename T>
void CenterMass<T>::operator ()(vector<Dvect> & v) {
	try {if(v.size() != cm.size()) throw " Center of mass dimensions differ !";}
	catch(const char * s) {cout << s <<endl; exit(1);}
    for(size_t n=0;n<v.size();n++){
      cm[n]=v[n];
    }
  };
template <typename T>
void CenterMass<T>::operator ()(vector<Dvect> & v, vector<Quaternion> & Q) {
	(*this)(v);
	try {if(Q.size() != Q_t.size()) throw " Quaternions dimensions differ !";}
	catch(const char * s) {cout << s <<endl; exit(1);}
	for(auto n=0;n<Q.size();n++){
		Q_t[n]=Q_t[n]*Q[n];

	}
}
template <typename T>
void CenterMass<T>::operator ()(vector<Dvect> & v,  vector<vector<Dvect> > & vm, vector<Quaternion> & Q,vector<Quaternion> & Qx ) {
	(*this)(v,Q);
	for(auto i=0;i<vm.size();i++){
		try{if(DoneOnce)
			if(cmm[i].size() != vm[i].size()) throw "The number of molecules in clusters changes from step to step!! ";}
		catch(const char * s){cout << s <<endl;exit(1);}
		cmm[i]=vm[i];
	}
	DoneOnce=true;
}

template class CenterMass<float>;
template class CenterMass<double>;
