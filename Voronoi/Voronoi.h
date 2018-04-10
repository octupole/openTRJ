/*
 * Voronoi.h
 *
 *  Created on: Jun 20, 2011
 *      Author: marchi
 */

#ifndef VORONOI_H_
#define VORONOI_H_
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include "Atoms.h"

#include <cstring>
#include <set>
#include "Array.h"
#include <fstream>
#include <iomanip>
#include "Topol.h"
#include "MyUtilClass.h"
#include "VoronoiSetter.h"
#include "Percolation.h"
#include "TrackOfstream.h"
#include "Split.h"
#include "voro++.hh"


using namespace Array;
using namespace MATRIX;
using namespace DVECT;
using Matrix=MATRIX::MMatrix<double>;
using namespace std;
using namespace voro;
using namespace Topol_NS;
using namespace std;
const int NNN=26;
const char h[]="Hxxxxx";
const int MAXSTR{7};
template <typename T>
inline char * as_byte(T & y){
	return reinterpret_cast<char *> (&y);
}
namespace Voro{
struct tArea{
	int n;
	double a;
	tArea(int m){n=m;a=0.0;}
	tArea(int m, double b){n=m;a=b;}
};

struct tAcomp{
	bool operator()(const tArea & x, const tArea & y) const {return x.n < y.n;};
};
class Voronoi {
protected:
	Matrix co{0};
	Matrix oc{0};
	static int nresid,nr,nc;
	int nWaters{0};
	static ios::streampos sizeHeader,sizeBody;
	bool readBinary{false};
	bool writeBinary{false};
	vector<int> SelectedResidues;
	static float time;
	vector<int> types; // Type number per atom according to ResidueTypes::_Residue collection
	vector<int> typesResidueMask;
	vector<int> atTypes;// Type of residue per atom
	vector<string> TypesName; // List of the residue types
	vector<double> Vol;
	vector<double> Rdii;
	array1<double> Vols;
	array2<double> area;
	vector<double> VolClusters,AreaClusters;
	vector<vector<double>> SurfaceClusters;
	container_periodic_poly * Mycon{nullptr};
	particle_order * porder{nullptr};
	static vector<string> label;
	vector<int> cindex;
	vector<vector<int> > CIndex;
	vector<vector<int> > Neighs;
	vector<vector<int> > NeighsRd;
	vector<vector<int> > * Neighsx;

	vector<vector<double> > Surface;
	vector<int> RealResidue;
	vector<string> Residue;
	vector<vector<int>> wShells;
	vector<vector<int>> Clusters;
	vector<int> atClusters;
	double VolCell{0};
	void gather(vector<int> & y);
	virtual void WriteIt(std::ofstream &)=0;
	virtual void ReadIt(std::ifstream &)=0;
	virtual void __compShell()=0;
	virtual void __searchNeighs(int a,int b)=0;
	virtual void __computeAggregate()=0;
	void dmpString(ofstream &,string &);
	void rdString(ifstream &,string &);
	template <typename T>
	void dmpVector(ofstream & f, vector<T> & v);
	template <typename T>
	void dmpVector(ofstream & f, vector<vector<T>> & v);
	void dmpVector(ofstream & f, vector<string> & v);
	template <typename T>
	void rdVector(ifstream & , vector<T> &);
	template <typename T>
	void rdVector(ifstream & , vector<vector<T>> &);
	void rdVector(ifstream & , vector<string> &);
	void bPrintBody(ofstream &);

	void bPrintHeader(ofstream &);
	void bReadHeader(ifstream &);
	void bReadBody(ifstream &);

	Voronoi();
public:
	Voronoi(Topol &,bool);
	Voronoi(Topol &,bool, Parallel::NewMPI *);
	Voronoi(ifstream & x);
	virtual void WriteLastJSON(std::ofstream &)=0;
	template <typename T>
	void Start(float, Atoms<T> &);
	void doVoro__();
	virtual void getData();

	const Matrix & getCO(){return co;}
	static float gTime(){return time;}
	size_t nClusters(){return Clusters.size();}
	void testVol();
	virtual ~Voronoi(){
		delete Mycon;
		delete porder;
	};
	string & name(int i) const {return label[i];}
	void setTypes(int,const AtomIndex [],const string *);
	double & getVolR(int n){return Vols[n];};
	int VoronoiCellID(const double & xc,const double & yc,const double & zc){
		double rx,ry,rz;
		int iD;
		try{if(!Mycon->find_voronoi_cell(xc,yc,zc,rx,ry,rz,iD))
			throw "Error while finding Voronoi cell ";}
		catch(const char * s) {std::cout << s << std::endl;exit(1);}
		return cindex[iD];
	}
	int VoronoiCellID(const double & xc,const double & yc,const double & zc,
			double & rx, double & ry, double & rz){
		int iD;
		try{if(!Mycon->find_voronoi_cell(xc,yc,zc,rx,ry,rz,iD))
			throw "Error while finding Voronoi cell ";}
		catch(const char * s) {std::cout << s << std::endl;exit(1);}
		return cindex[iD];
	}
	int getTypes(int i) {return atTypes[i];};
	int getTypesRes(int nc) {
		int First=CIndex[nc][0];
		return getTypes(First);
	}
	string & getTypesName(int i) {return TypesName[types[i]];};
	friend std::ofstream & operator<<(std::ofstream & fout, Voronoi  & y){
		y.WriteIt(fout);
		return fout;
	}
	friend std::ifstream & operator>>(std::ifstream & fin, Voronoi  & y){
		y.ReadIt(fin);
		return fin;
	}
};
}
#endif /* VORONOI_H_ */
