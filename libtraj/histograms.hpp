/*
 * histograms.hpp
 *
 *  Created on: Aug 16, 2013
 *      Author: marchi
 */

#ifndef HISTOGRAMS_HPP_
#define HISTOGRAMS_HPP_
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "Ftypedefs.h"



using namespace std;

const string MyAxis[DIM]={"X","Y","Z"};

struct VolumeData{
	int nX,nY,nZ;
	double Xmax,Xmin,Xstride;
	double Ymax,Ymin,Ystride;
	double Zmax,Zmin,Zstride;
	vector<double> data;
	VolumeData(int nx, int ny, int nz, double xstride,double ystride,double zstride ): nX(nx),nY(ny),nZ(nz),Xstride(xstride)
		,Ystride(ystride),Zstride(zstride), Xmin(0.0), Ymin(0.0), Zmin(0.0){
		Xmax=static_cast<double>(nX-1)*Xstride;
		Ymax=static_cast<double>(nY-1)*Ystride;
		Zmax=static_cast<double>(nZ-1)*Zstride;
		data=vector<double>(nX*nY*nZ,0.0);
	}
	VolumeData(int nx, int ny, int nz): nX(nx),nY(ny),nZ(nz),Xstride(0.0)
		,Ystride(0.0),Zstride(0.0), Xmin(0.0), Ymin(0.0), Zmin(0.0), Xmax(0.0), Ymax(0.0), Zmax(0.0){
		data=vector<double>(nX*nY*nZ,0.0);
	}
	VolumeData(): nX(0),nY(0),nZ(0),Xstride(0.0)
		,Ystride(0.0),Zstride(0.0), Xmin(0.0), Ymin(0.0), Zmin(0.0), Xmax(0.0), Ymax(0.0), Zmax(0.0){
	}
	VolumeData & operator+=(VolumeData y){
		try{
			if(nX != 0 && nY!=0 && nZ!=0)
				if(nX != y.nX || nY != y.nY || nZ != y.nZ) throw string("Can't accumulate VolumeData, "
						"dimensions are not identical.");
			else{
				nX=y.nX;nY=y.nY;nZ=y.nZ;
				data=vector<double>(nX*nY*nZ,0.0);
			}
		} catch(const string & s){
			cout << s <<endl;
			exit(1);
		}
		Xmax=y.Xmax;Ymax=y.Ymax;Zmax=y.Zmax;
		Xmin=y.Xmin;Ymin=y.Ymin;Zmin=y.Zmin;
		Xstride=y.Xstride;Ystride=y.Ystride;Zstride=y.Zmax;
		for(size_t o=0; o < data.size();o++)
			data[o]+=y.data[o];

		return *this;
	}
	friend ostream & operator<<(ostream & fout, VolumeData y){
		fout << fixed << "[ ";
		fout << fixed << setw(8) << setprecision(4) << y.Xmin <<":";
		fout << fixed << setw(8) << setprecision(4) << y.Xstride <<":";
		fout << fixed << setw(8) << setprecision(4) << y.Xmax <<" ] ";
		fout << fixed << "[ ";
		fout << fixed << setw(8) << setprecision(4) << y.Ymin <<":";
		fout << fixed << setw(8) << setprecision(4) << y.Ystride <<":";
		fout << fixed << setw(8) << setprecision(4) << y.Ymax <<" ] ";
		fout << fixed << "[ ";
		fout << fixed << setw(8) << setprecision(4) << y.Zmin <<":";
		fout << fixed << setw(8) << setprecision(4) << y.Zstride <<":";
		fout << fixed << setw(8) << setprecision(4) << y.Zmax <<" ] ";
		fout <<endl;
		fout << fixed << setw(8) << y.nX << " " << y.nY << " " << y.nZ << " ";
		fout <<endl;
		int na=0;
		for(int o=0;o< y.nX;o++){
			for(int p=0;p < y.nY;p++){
				for(int q=0;q < y.nZ;q++){
					fout << fixed << setw(14) << right << scientific << setprecision(4) << y.data[na] << " ";
					na++;
				}
			}
		}
		fout << endl;
		return fout;
	}
};
struct hist1D{
	double hisX;
	double nhis;
	hist1D():hisX(0.0),nhis(0){};
	hist1D(double his0,double nhis0): hisX(his0),nhis(nhis0){};
	hist1D & operator+=(hist1D a){
		hisX+=a.hisX; nhis+=a.nhis;
		return *this;
	}
	hist1D & operator/=(double a){
		hisX/=a;
		return *this;
	}
	hist1D & operator*=(double a){
		hisX*=a;
		return *this;
	}
	bool isZero(){
		return nhis == 0;
	}
	double Ratio(){
		if(!isZero()) return hisX/static_cast<double> (nhis);
		else return 0.0;
	}

	hist1D & operator=(const hist1D & y){
		hisX=y.hisX;
		nhis=y.nhis;
		return *this;
	}
	double gHisx(){return hisX;}
};
struct Histogram1D{
private:
	virtual void WriteIt(ostream & );
public:
	vector<hist1D> hist;
	double dx{0.0},cut{0.0};
	int HisX{0};
	double Factor{1.0/unit_nm};
	string Label{" "};
	Histogram1D(){};
	Histogram1D(double,double);
	void operator()(double,double);

	virtual void setLabel(string s){Label=s;}
	void setFactor(double fact){Factor=fact;}

	virtual hist1D & operator[](size_t);
	Histogram1D  operator*(double );
	Histogram1D  operator*(size_t );
	Histogram1D & operator=(const Histogram1D &);
	bool operator==(const Histogram1D &);

	virtual Histogram1D & operator+=(const Histogram1D &);
	virtual void clear();
	friend ostream & operator<<(ostream &, Histogram1D & );

	virtual ~Histogram1D(){};
};

struct Histogram2D: public Histogram1D {
	vector<Histogram1D> His2D;
	int HisY;
	double dy;
	double cutz;
	static string Tag[DIM];
	Histogram2D(): Histogram1D(), dy(dx), cutz(0.0), HisY(0){};
	Histogram2D(double dx0,double cutxy0, double cutz0 ):
		Histogram1D(dx0,cutxy0), dy(dx0), cutz(cutz0), HisY(0){
		int nh=static_cast<int>(cut/dy);
		HisY=nh;
		His2D=vector<Histogram1D>(2*nh+1,dynamic_cast<Histogram1D &>(*this));
	};
	Histogram2D & operator+=(const Histogram2D & z){
		dx=z.dx;
		cut=z.cut;
		Label=z.Label;
		HisX=z.HisX;
		HisY=z.HisY;
		HisX=z.HisX;
		cutz=z.cutz;
		dy=z.dy;
		try{
			if(His2D.size() == 0){
				His2D=vector<Histogram1D>(z.His2D.size(),Histogram1D());
			} else if(His2D.size() != z.His2D.size()) throw string("The His2D dimension are not identical ");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
		for(size_t o=0;o<His2D.size();o++){
			His2D[o]+=z.His2D[o];
		}
		return *this;
	}
	Histogram2D & operator=(const Histogram2D & y){
		His2D=y.His2D;
		HisY=y.HisY;
		HisX=y.HisX;
		cutz=y.cutz;
		dy=y.dy;
		this->Histogram1D::operator=(*dynamic_cast<const Histogram1D *>(&y));
		return *this;
	}
	void operator()(double dx0, double cutxy0, double cutz0 ){
		dy=dx0;cutz=cutz0;
		int nh=static_cast<int>(cut/dy);
		dynamic_cast<Histogram1D &>(*this)(dx,cutxy0);
		HisY=nh;
		His2D=vector<Histogram1D>(2*nh+1,dynamic_cast<Histogram1D &>(*this));
	}
	Histogram1D & operator[](int n){
		return His2D[HisY+n];
	}
	size_t size(){return His2D.size();}
	void static sTags(vector<int> y){
		for(int o=0;o<DIM;o++){
			Tag[o]=MyAxis[y[o]];
		}
	}
	friend ostream & operator<<(ostream & fout, Histogram2D y){
		fout << fixed << "[ ";
		fout << fixed << setw(8) << setprecision(2) << -static_cast<double>(y.HisX)*y.dx/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << y.dx/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << static_cast<double>(y.HisX)*y.dx/unit_nm <<" ] ";
		fout << fixed << "[ ";
		fout << fixed << setw(8) << setprecision(2) << -static_cast<double>(y.HisY)*y.dy/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << y.dy/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << static_cast<double>(y.HisY)*y.dy/unit_nm <<" ] ";
		fout <<endl;
		fout << fixed << setw(8) << 2*y.HisX +1 << " " << 2*y.HisY+1 << " ";
		fout <<endl;
		for(int o=-y.HisX;o<= y.HisX;o++){
			double ddx=y.dx*static_cast<double>(o);
			for(int p=-y.HisY;p<=y.HisY;p++){
				double ddy=y.dy*static_cast<double>(p);
				double f=y[o][p].Ratio();
				fout << fixed << setw(14) << right << scientific << setprecision(7) << f << " ";
			}
		}
		fout << endl;
		return fout;

	}
	virtual ~Histogram2D(){};
};
struct Histogram3D: public Histogram2D {
	vector<Histogram2D> His3D;
	int HisZ;
	Histogram3D(): Histogram2D(), HisZ(0){};
	Histogram3D(double dx0,double cut0):
		Histogram2D(dx0,cut0,cut0), HisZ(0){
		int nh=static_cast<int>(cut/dx);
		HisZ=nh;
		His3D=vector<Histogram2D>(2*nh+1,dynamic_cast<Histogram2D &>(*this));
	};

	Histogram3D & operator=(const Histogram3D & y){
		His3D=y.His3D;
		HisZ=y.HisZ;
		HisY=y.HisY;
		HisX=y.HisX;
		cut=y.cut;
		dy=y.dy;
		dx=y.dx;
		this->Histogram2D::operator=(*dynamic_cast<const Histogram2D *>(&y));
		return *this;
	}
	void operator()(double dx0, double cut0 ){
		dy=dx0;cut=cut0;
		int nh=static_cast<int>(cut/dx);
		dynamic_cast<Histogram2D &>(*this)(dx,cut0,cut0);
		HisZ=nh;
		His3D=vector<Histogram2D>(2*nh+1,dynamic_cast<Histogram2D &>(*this));
	}
	Histogram2D & operator[](int n){
		return His3D[HisZ+n];
	}
	Histogram3D & operator+=(const Histogram3D & z){
		HisZ=z.HisZ;
		HisY=z.HisY;
		HisX=z.HisX;
		cut=z.cut;
		dy=z.dy;
		dx=z.dx;
		try{
			if(His3D.size() == 0){
				His3D=vector<Histogram2D>(z.His3D.size(),Histogram2D());
			} else if(His3D.size() != z.His3D.size()) throw string("The His2D dimension are not identical ");
		}catch(const string & s){
			cout << s << endl;
			exit(1);
		}
		for(size_t o=0;o<His3D.size();o++){
			His3D[o]+=z.His3D[o];
		}
		return *this;
	}
	size_t size(){return His3D.size();}

	friend ostream & operator<<(ostream & fout, Histogram3D y){
		fout << fixed << "[ ";
		fout << fixed << setw(8) << setprecision(2) << -static_cast<double>(y.HisX)*y.dx/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << y.dx/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << static_cast<double>(y.HisX)*y.dx/unit_nm <<" ] ";
		fout << fixed << "[ ";
		fout << fixed << setw(8) << setprecision(2) << -static_cast<double>(y.HisY)*y.dx/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << y.dy/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << static_cast<double>(y.HisY)*y.dx/unit_nm <<" ] ";
		fout << fixed << "[ ";
		fout << fixed << setw(8) << setprecision(2) << -static_cast<double>(y.HisZ)*y.dx/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << y.dy/unit_nm <<":";
		fout << fixed << setw(8) << setprecision(2) << static_cast<double>(y.HisZ)*y.dx/unit_nm <<" ] ";
		fout <<endl;
		fout << fixed << setw(8) << 2*y.HisX +1 << " " << 2*y.HisY+1 << " " << 2*y.HisZ+1 << " ";
		fout <<endl;
		for(int o=-y.HisX;o<= y.HisX;o++){
			for(int p=-y.HisY;p<=y.HisY;p++){
				for(int q=-y.HisZ;q<=y.HisZ;q++){
					double f=y[o][p][q].Ratio();
					fout << fixed << setw(14) << right << scientific << setprecision(7) << f << " ";
				}
			}
		}
		fout << endl;
		return fout;

	}
	virtual ~Histogram3D(){};
};

#endif /* HISTOGRAMS_HPP_ */
