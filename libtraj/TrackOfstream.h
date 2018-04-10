/*
 * TrackOfstream.h
 *
 *  Created on: Nov 22, 2017
 *      Author: marchi
 */

#ifndef LIBTRAJ_TRACKOFSTREAM_H_
#define LIBTRAJ_TRACKOFSTREAM_H_
#include <ostream>
#include <fstream>
#include <string>

using std::string;
using std::ofstream;
using std::ios;
class TrackOfstream {
	ofstream * out{nullptr};
	ios::off_type countW{0};
public:
	TrackOfstream(ofstream * x): out{x}{};
	void write(const char * c,ios::off_type size){
		out->write(c,size);
		countW+=size;
	}
	void write(char * c,ios::off_type size){
		out->write(c,size);
		countW+=size;
	}
	ios::off_type gCountW(){return countW;}
	void set0CountW(){countW=0;}

	virtual ~TrackOfstream(){
		if(out) delete out;
	};
};
/*
 #include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <sstream>
using namespace std;
template <typename T>
inline char * as_byte(T & y){
	return reinterpret_cast<char *> (&y);
}

class streams{
	ofstream * out;
	ios::off_type countW{0};
public:
	streams(ofstream * x){out=x;};
	void write(const char * c,ios::off_type size){
		out->write(c,size);
		countW+=size;
	}
	void write(char * c,ios::off_type size){
		out->write(c,size);
		countW+=size;
	}
	ios::off_type gCountW(){return countW;}
};

const size_t NMAX{1234},NCYCLE{10};
int main(){
	 // get size of file

	ofstream fout;
	streams Fout{&fout};
	ifstream fin;
	stringstream iss;

	string myFile1{"guga1.bin"};
	vector<double> v(NMAX,0.0);
	vector<double> vv(NMAX,-1.0);
	iss.write(as_byte(v[0]),sizeof(v[0])*NMAX);
	size_t size=sizeof(v[0])*NMAX;
	cout << size<<endl;
	fout.open(myFile1,ios::out|ios::binary);

	Fout.write(iss.str().c_str(),sizeof(v[0])*NMAX);
	for(size_t o{0};o< NCYCLE;o++)
		fout.write(iss.str().c_str(),sizeof(v[0])*NMAX);

	fout.seekp(2*size,fout.beg);
	fout.write(as_byte(vv[0]),sizeof(v[0])*NMAX);
	fout.close();
	fin.open(myFile1,ios::in|ios::binary);
	vector<double> vg(NMAX);
	cout << NMAX<<endl;
	fin.seekg(2*size-sizeof(vg[0]),fin.beg);
	fin.read(as_byte(vg[0]),sizeof(vg[0]));
	cout << " " << vg[0]<<endl;
	fin.read(as_byte(vg[0]),sizeof(vg[0]));
	cout << " " << vg[0]<<endl;
	cout << Fout.gCountW()<<endl;



}


 */
#endif /* LIBTRAJ_TRACKOFSTREAM_H_ */
