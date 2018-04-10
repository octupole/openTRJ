/*
 * ExecuteVoronoi.cpp
 *
 *  Created on: Nov 16, 2017
 *      Author: marchi
 */

#include "../Execute/ExecuteVoronoi.h"

#include "Atoms.h"
#include "Voronoi.h"
#include "VoronoiMicellesJSON.h"
#include "VoronoiBinary.h"

namespace Voro {

template <typename T>
size_t ExecuteVoronoi<T>::nnx=1;
template <typename T>
size_t ExecuteVoronoi<T>::nny=1;
template <typename T>
size_t ExecuteVoronoi<T>::nnz=1;
template <typename T>
Parallel::NewMPI * ExecuteVoronoi<T>::CurrMPI{nullptr};
template <typename T>
void ExecuteVoronoi<T>::__Print(ostream & y){}

template <typename T>
ExecuteVoronoi<T>::ExecuteVoronoi(trj::TrjRead & MyIn) {
	__SetUp(MyIn);

	bOnce=MyIn.bbOnce();
	nnx=MyIn.gnnx();
	nny=MyIn.gnny();
	nnz=MyIn.gnnz();
	nstart=MyIn.gnstart();
	nend=MyIn.gnend();
	nskip=MyIn.gnskip();
	bTest=MyIn.bbTestVol();

	if(finx){
		size_t TotFrame=finx->gFrameStep();
		try{
			if(int(TotFrame) <nend) {
				stringstream ss0,ss1;
				ss0<< nend; ss1<< TotFrame;
				throw string("\n\nThe requested end-point \"")+ss0.str()
									+string("\" goes beyond the last frame of the \ntrajectory \""+ss1.str()+string("\"!!\n\n"));
			}
		}catch(const string & s){
			cout << s << endl;
			Finale::Finalize::Final();
		}
	}
	try{
		if(finx && nend-nstart+1 < nskip) throw string("\nNumber of selected steps is smaller than -skip parameter. Change and rerun.\n");
	}catch(const string & s) {cout << s <<endl;Finale::Finalize::Final();}

	if(fin1x){
		cout << JSONOutput<<endl;
		if(JSONOutput)
			vor=new VoronoiBinary<VoronoiMicellesJSON>(*fin1x);
		else
			vor=new VoronoiBinary<VoronoiMicelles>(*fin1x);

		ofstream & fout=*foutx;

		CurrMPI->Barrier();
		Comms=new Parallel::FComms(CurrMPI,fout,fileout,nstart,nend,1e8,1);
		nstart=Comms->getStart();
		nend=Comms->getEnd();
	}
}
template <typename T>
ExecuteVoronoi<T>::ExecuteVoronoi(trj::TrjRead & MyIn, Topol & Topology):
 	 ExecuteVoronoi<T>::ExecuteVoronoi(MyIn){
	Top=&Topology;
	ios::streampos len;
	HeaderTrj header;
// Read header of dcd file
	try{
		if(finx) {
			finx->seekg(0,"end");
			len=finx->tellg();
			finx->seekg(0,"beg");
			*finx>>header;
			try{
				if(!header.check(Top->Size())) {
					stringstream ss0,ss1;
					ss0<< Top->Size();
					ss1<<header.getNatoms();
					throw string("Number of atoms in the pdb ("+ss0.str()+") and trajectory files ("
							+ss1.str()+") does not match!");}
			}
			catch(const string & s){cout << s<<endl;Finale::Finalize::Final();}
			ofstream & fout=*foutx;

			CurrMPI->Barrier();

			Comms=new Parallel::FComms(CurrMPI,fout,fileout,nstart,nend,header.getNFR(),nskip);
			nstart=Comms->getStart();
			nend=Comms->getEnd();
		}else{
			ofstream & fout=*foutx;
			CurrMPI->Barrier();
			nstart=1;nend=1;nskip=1;
			if(CurrMPI->Get_Size() > 1) throw string("Cannot compute Voronoi from a PDB file in Parallel!");
			Comms=new Parallel::FComms(CurrMPI,fout,fileout,nstart,nend,header.getNFR(),nskip);
			nstart=Comms->getStart();
			nend=Comms->getEnd();
		}
	} catch(const string & s){cout << s<<endl;
	Finale::Finalize::Final();exit(1);}

	if(JSONOutput){
		if(binOutput){
			vor=new VoronoiBinary<VoronoiMicellesJSON>(*foutx,Topology,bHyd, CurrMPI);
		}else{
			vor=new VoronoiMicellesJSON(Topology,bHyd);
		}
	}else{
		if(binOutput){
			vor=new VoronoiBinary<VoronoiMicelles>(*foutx,Topology,bHyd, CurrMPI);
		} else{
			vor=new VoronoiMicelles(Topology,bHyd);
		}
	}
	Percolation<T>::setPercoCutoff(MyIn.gPercoCutoff());

	Clustering=MyIn.bbClust();
	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();

};
template <typename T>
void ExecuteVoronoi<T>::operator()(Atoms<T> * atx){

	if(!atx) {
		__RunPost();
		return;
	}
	if(finx)
		__RunTrajectory(atx);
	else
		__RunPDB(atx);
}

template <typename T>
void ExecuteVoronoi<T>::__RunPost(){
	myiterators::IteratorVoronoi iter_vor(vor,fin1x,nstart,nend);
	while((++iter_vor).isReferenced()){
		stringstream ss;
		Voronoi * vorA=iter_vor();
		vorA->getData();
		Comms->getStream() << *vorA;
		if(bTest) vor->testVol();
		switch(vor->nClusters()){
		case 0:
			break;
		case 1:
			ss << "    " <<fixed << setw(4) << vor->nClusters()<<" cluster  <-----";
			break;
		default:
			ss<< "    " << fixed << setw(4) << vor->nClusters()<<" clusters <-----";
		}



		cout << fixed << setw(5) << "----> Time Step " << vor->gTime() << ss.str()<<"\n";
	}
	vor->WriteLastJSON(Comms->getStream());
	Comms->appendStreams();
	if(bDel) Comms->removeFiles();
	CurrMPI->Barrier();
	CurrMPI->~NewMPI();
	cout << "\nProgram completed: Output data written to " + fileout << "\n\n";
}
template <typename T>
void ExecuteVoronoi<T>::__RunTrajectory(Atoms<T> * atmx){

	myiterators::IteratorAtoms<T> iter_atm(atmx,finx,nstart,nend,nskip);

	while((++iter_atm).isReferenced()){
		stringstream ss;
		Atoms<T> * atmA=iter_atm();

		float ntime=atmA->getTime();
		int nClusters{0};
		atmA->setTopol(*Top);
		if(Clustering){
			static struct Once{
				Once(Atoms<T> * atmA, Topol_NS::Topol * myTop, bool JSON){
					if(JSON)
							atmA->template SetupPercolate<Enums::JSON>(*myTop);
						else
							atmA->template SetupPercolate<Enums::noJSON>(*myTop);
				}
			} _Once(atmA, Top, this->JSONOutput);
			if(bOnce){
				static struct Once_p{int nClusters{0};Once_p(Atoms<T> *atmA){nClusters=atmA->Percolate();}} __Once_p(atmA);
				nClusters=__Once_p.nClusters;
			}else {
				nClusters=atmA->Percolate();
			}
		}
		switch(nClusters){
		case 0:
			break;
		case 1:
			ss << "    " <<fixed << setw(4) << nClusters<<" cluster  <-----";
			break;
		default:
			ss<< "    " << fixed << setw(4) << nClusters<<" clusters <-----";
		}

		vor->Start(ntime,*atmA);
		vor->doVoro__();
		vor->getData();

		Comms->getStream() << *vor;

		if(bTest) vor->testVol();
		cout << fixed << setw(5) << "----> Time Step " << ntime << ss.str()<<"\n";
	}

	vor->WriteLastJSON(Comms->getStream());
	Comms->appendStreams();
	if(bDel) Comms->removeFiles();
	Comms->closeStream();   // close stream!!!
	CurrMPI->Barrier();
	CurrMPI->~NewMPI();
	cout << "\nProgram completed: Output data written to " + fileout << "\n\n";

}
template <typename T>
void ExecuteVoronoi<T>::__RunPDB(Atoms<T> * atm){
	bool bTestVol{true};
	vector<string> data;
	fpdb->clear();
	fpdb->seekg(0);
	for(string str;getline(*fpdb,str);){
		data.push_back(str);
	}
	atm->pdb(data);
	float ntime=atm->getTime();
	vor-> Start(ntime,*atm);
	vor->doVoro__();
	vor->getData();
	Comms->getStream() << *vor;

	if(bTestVol) vor->testVol();
	Comms->appendStreams();
	if(bDel) Comms->removeFiles();
	CurrMPI->Barrier();
	CurrMPI->~NewMPI();
	cout << "\nProgram completed: Output data written to " + fileout << "\n\n";
}


template <typename T>
void ExecuteVoronoi<T>::__SetUp(trj::TrjRead & MyIn){
	fileout=MyIn.gfileout();

	bDel=MyIn.bbDel();
	bHyd=MyIn.bbHyd();
	binOutput=MyIn.bbOutBin();
	JSONOutput=MyIn.bbOutJSON();
	fileout_bin=MyIn.gfileout_bin();

	fpdb=MyIn.gFpdb();
	finx=MyIn.gFinx();
	fin1x=MyIn.gFin1();
	foutx=MyIn.gFoutx();

	fidb=MyIn.gFidb();
	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
	foutx=new ofstream();


	/*
	 * Define and dimension RhoSaxs and Saxs classes
	 */

	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
	Clustering=MyIn.bbClust();
}

template <typename T>
ExecuteVoronoi<T>::~ExecuteVoronoi() {}

template class ExecuteVoronoi<float>;
template class ExecuteVoronoi<double>;

} /* namespace Voronoi */
