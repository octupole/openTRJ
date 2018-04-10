/*
 * LCells.cpp
 *
 *  Created on: Dec 18, 2014
 *      Author: marchi
 */

#include "LCells.h"

template <typename T>
LCells<T>::LCells(): nr(0){
	// TODO Auto-generated constructor stub

}
template<typename T>
bool LCells<T>::test(){
	if(nc[XX] < 5) return false;
	if(nc[YY] < 5) return false;
	if(nc[ZZ] < 5) return false;
	return true;
}

template <typename T>
void LCells<T>::Init(Matrix & co0, const vector<Dvect> & y){
	T Rcut0=Rcut;
	x=y;
	co=co0;
	oc=co.Inversion();
	nr=x.size();
	try{
		if(co[XX][XX] < Rcut*2.0) throw string{"Box is too small to run with neighbor lists."};
	}catch(const string & s){
		cout << s<<endl;
		exit(1);
	}
	if(co[XX][XX] > Rmax){
		if(nc[XX] < 0){
			nc[XX]=static_cast<int>(co[XX][XX]/Rcut0);
			nc[YY]=static_cast<int>(co[YY][YY]/Rcut0);
			nc[ZZ]=static_cast<int>(co[ZZ][ZZ]/Rcut0);
		}
		if(Chainp.size()) Chainp.clear();
		Chainp=vector<vector<vector<vectint> > >(nc[XX]);
		for(int m=0;m<nc[XX];m++){
			Chainp[m]=vector<vector<vectint> >(nc[YY]);
			for(int n=0;n<nc[YY];n++)
				Chainp[m][n]=vector<vectint>(nc[ZZ]);
		}
	}
	nnl.clear();
	nnl=vector<vectint>(nr);
}

template <typename T>
void LCells<T>::Index(){
	if(nc[XX] < 0 && nc[YY] < 0 && nc[ZZ] < 0) return;
	T Rcut0=Rcut;
	T sqcut=Rcut0*Rcut0;
	vector<int> nn=vector<int>(DIM,0);
	if(indx.size()) indx.clear();

	indx.push_back(nn);
	int imax=0,jmax=0,kmax=0;
	for(int i=-nc[XX];i<nc[XX];i++){
		nn[XX]=i;
		for(int j=-nc[YY];j<nc[YY];j++){
			nn[YY]=j;
			for(int k=-nc[ZZ];k<nc[ZZ];k++){
				nn[ZZ]=k;
				T dmin=Dist_ijk(i,j,k);
				bool v= (i == 0 && j == 0 && k == 0);
				if(dmin < sqcut && !v){
					indx.push_back(nn);
					imax=(imax < abs(i))?abs(i):imax;
					jmax=(jmax < abs(j))?abs(j):jmax;
					kmax=(kmax < abs(k))?abs(k):kmax;
				}
			}
		}
	}

	try{
		if( imax > (nn[XX])/2 || jmax > (nn[YY])/2 || kmax > (nn[ZZ])/2)
			throw " Increase Linked Cells dimensions ";
	}
	catch(const char * s){
		cout << s << endl;
		exit(0);
	}

}

template <typename T>
vector<vector<int> > & LCells<T>::List(bool UpperHalf){

	if(nc[XX] < 0 && nc[YY] < 0 && nc[ZZ] < 0) {
		for(size_t n=0;n< x.size();n++){
			Dvect X=co*x[n];
			for(auto o=UpperHalf?n+1:0;o<x.size();o++){
				if(n == o) continue;
				Dvect Y=co*x[o];
				T dist=Y.Dist(X,co,oc);
				if(dist<Rcut) {
					nnl[n].push_back(o);
				}

			}
		}
		return nnl;
	}

	T d[3]={static_cast<T>(nc[XX]),static_cast<T>(nc[YY]),static_cast<T>(nc[ZZ])};
	if(Cellp.size()){
		Cellp.clear();
	}
	Cellp=vector<vectint>(x.size());
	for(size_t i=0;i< x.size();i++){
		Cellp[i]=vector<int>(DIM,0);
		T x1 = x[i][XX];
		T y1 = x[i][YY];
		T z1 = x[i][ZZ];
		int mx = static_cast<int> (d[XX] * (x1 - rint(x1 - HALF)));
		int my = static_cast<int> (d[YY] * (y1 - rint(y1 - HALF)));
		int mz = static_cast<int> (d[ZZ] * (z1 - rint(z1 - HALF)));
		Cellp[i][XX]=(mx%nc[XX]+nc[XX])%nc[XX];
		Cellp[i][YY]=(my%nc[YY]+nc[YY])%nc[YY];
		Cellp[i][ZZ]=(mz%nc[ZZ]+nc[ZZ])%nc[ZZ];
		Chainp[mx][my][mz].push_back(i);
	}

	for(size_t n=0;n< x.size();n++){
		Dvect X=co*x[n];
		vectint i=Cellp[n];
		for(size_t p=0;p<indx.size();p++){
			int nx=indx[p][XX]+i[XX];
			nx=(nx%nc[XX]+nc[XX])%nc[XX];
			int ny=indx[p][YY]+i[YY];
			ny=(ny%nc[YY]+nc[YY])%nc[YY];
			int nz=indx[p][ZZ]+i[ZZ];
			nz=(nz%nc[ZZ]+nc[ZZ])%nc[ZZ];
			size_t mm=Chainp[nx][ny][nz].size();
			for(size_t m0=0;m0<mm;m0++){
				size_t m=Chainp[nx][ny][nz][m0];
				if(m == n) continue;
				if(UpperHalf && m<n) continue;
				Dvect Y=co*x[m];
				T dist=Y.Dist(X,co,oc);
				if(dist < Rcut) {
					nnl[n].push_back(m);
				}

			}
		}
	}
	return nnl;
}
template <typename T>
T LCells<T>::Dist_ijk(int ni,int nj, int nk){
	T dmin=1.0e+8;
	vector<T> d=vector<T>(DIM);
	for(size_t n=0;n<DIM;n++){
		d[n]=1.0/static_cast<T>(nc[n]);
	}
	vector<vector<int> > nv;
	for(int o=0;o<2;o++)
		for(int p=0;p<2;p++)
			for(int q=0;q<2;q++){
				int n[]={o,p,q};
				nv.push_back(vector<int>(n,n+DIM));
			}
	for(size_t n=0;n<nv.size();n++)
		for(size_t m=0;m<nv.size();m++){
			T lx=static_cast<T>(ni+nv[m][XX]-nv[n][XX])*d[XX];
			T ly=static_cast<T>(nj+nv[m][YY]-nv[n][YY])*d[YY];
			T lz=static_cast<T>(nk+nv[m][ZZ]-nv[n][ZZ])*d[ZZ];
			T mx=co[XX][XX]*lx+co[XX][YY]*ly+co[XX][ZZ]*lz;
			T my=co[YY][XX]*lx+co[YY][YY]*ly+co[YY][ZZ]*lz;
			T mz=co[ZZ][XX]*lx+co[ZZ][YY]*ly+co[ZZ][ZZ]*lz;
			T d=mx*mx+my*my+mz*mz;
			if(d < dmin){
				dmin=d;
			}

		}
	try{
		if(fabs(co[XX][YY]) > small || fabs(co[XX][ZZ]) >small || fabs(co[YY][ZZ]) > small) throw
				" Cannot yet work with non orthogonal axis! ";
	}
	catch(const char * s){
		std::cout << s << std::endl;
		exit(1);
	}
	return dmin;
}
template <typename T>
LCells<T>::~LCells() {
	// TODO Auto-generated destructor stub
}
template class LCells<float>;
template class LCells<double>;
