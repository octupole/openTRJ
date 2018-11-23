/*
 * ExecbSaxsTraj.cpp
 *
 *  Created on: Dec 16, 2015
 *      Author: marchi
 */

#include "ExecbSaxsTraj.h"


enum Model {Consecutive, Chunks};
Model MyModel{Chunks};
auto min3=[](Matrix co){
  double myMin=1e10;
  for(auto q=0;q<DIM;q++){
    if(myMin > 5.0*co[q][q]) myMin=5.0*co[q][q];
  }
  return myMin;
};


ExecbSaxsTraj::ExecbSaxsTraj(trj::TrjRead & MyIn, Topol & Topology):ExecbSaxs(MyIn), Top(&Topology){
	bool bFluct=MyIn.bbFluct();
	if(!bFluct) MyModel=Consecutive;
	bOnce=MyIn.bbOnce();
	nnx=MyIn.gnnx();
	nny=MyIn.gnny();
	nnz=MyIn.gnnz();
	nstart=MyIn.gnstart();
	nend=MyIn.gnend();
	nskip=MyIn.gnskip();
	try{
		if(nnx ==1 || nny == 1|| nnz == 1) throw string("Grid dimensions are not set!");
	}catch(const string & s){
		cout << s <<endl;
		Finale::Finalize::Final();
	}
	if(finx){
		size_t TotFrame=finx->gFrameStep();
		try{
			if(TotFrame <nend) {
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

	if(CurrMPI->AmI_Parallel()){

		auto MyTaskNo=CurrMPI->Get_Rank();
		try{
			if(finx && (nend-nstart+1)/nskip < CurrMPI->Get_Size()) throw string(" The number of CPUs is larger than the number "
					"of the trajectory steps. Change and rerun. ");
		}catch(const string & s) {cout << s <<endl;Finale::Finalize::Final();}
		if(MyModel == Consecutive){
			auto nstart_init=nstart;
			auto ntask=CurrMPI->Get_Size();
			auto ntot=(nend-nstart+1);
			nstart=nstart_init+MyTaskNo*ntot/ntask;
			nend=nstart+ntot/ntask;
		}
		else if(MyModel == Chunks){
			auto ntot=(nend-nstart+1)/nskip;
			auto ntask=CurrMPI->Get_Size();
			auto nChunk=(ntot/ntask)*nskip;
			nend=nstart+nChunk*(MyTaskNo+1)-1;
			nstart=nstart+nChunk*MyTaskNo;
		}

	}
	__SetUp(MyIn);

};
void ExecbSaxsTraj::__MakeDummyAtoms(){
	if(Patm) return;
	TopolPDB topPDB;
	Topol MyTop;
	vector<string> data;
	string str{"ATOM      1  DU  DUM     1       0.000   0.000   0.000  1.00  0.00            "};
	data.push_back(str);


	topPDB(data);
	MyTop(topPDB,false);
	TopPBC=new Topol(MyTop);

	int natoms=MyTop.Size();
	vector<string> SelRes{"DUM"};
	vector<string> Ref{"DUM"};

	Patm=new MAtoms(diffk);

	Patm->setDim(natoms);
	Patm->pdb(data);
	Patm->initLists(topPDB, SelRes);
	Patm->InitSelection<Selection>(SelRes,MyTop);
	Patm->InitSelection<Reference>(Ref,MyTop);

}

void ExecbSaxsTraj::operator()(MAtoms * atm){

	if(finx)
		__RunTrajectory(atm);
	else
		__RunPDB(atm);
}

void ExecbSaxsTraj::__SetUp(trj::TrjRead & MyIn){
	ExecbSaxs::__SetUp(MyIn);
	SuperCell0=MyIn.gSuperCellSide();

	MyOrder=WhatToDo==Enums::ELDENS?4:MyIn.gMyOrder();
	bSaxsBSP=WhatToDo==Enums::ELDENS?false:MyIn.bbSaxsBSP();
	bFluct=MyIn.bbFluct();
	bFixed=MyIn.bbFixed();
	nacc=MyIn.gnacc();
	bDebye=MyIn.bbDebye();
	bDirect=MyIn.bbDirect();
	exPadding=MyIn.WhichPadding;


	if(exPadding == Enums::Periodic) Filter.Allocate(nnx,nny,nnz);

	ios::streampos len;
	HeaderTrj header;
// Read header of dcd file
	if(finx) {
		finx->seekg(0,"end");
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
	}
	/*
	 * Define and dimension RhoSaxs and Saxs classes
	 */
	if(!bDebye && !bDirect){
		if((bSaxsBSP || !bFluct ) && WhatToDo!=Enums::ELDENS){
			Rho_ex=new RhoSaxsBSP;
			if(!bFluct && !bFixed)
				MySaxs=new SaxsBSPStatic(nacc,MyOrder,Myd,Mycut);
			else if(bFixed)
				MySaxs=new SaxsBSPfixed(MyOrder,Myd,Mycut);
			else
				MySaxs=new SaxsBSP(MyOrder,Myd,Mycut);
		} else {
			if(MyOrder>1) Rho_ex=new RhoSaxsLI;
			else Rho_ex=new RhoSaxs;
			MySaxs=new Saxs(MyOrder,Myd,Mycut);
		}
	}else if(bDebye){
		MySaxs=new SaxsDebye(MyOrder,Myd,Mycut);
	} else if(bDirect){
		MySaxs=new SaxsDir(MyOrder,Myd,Mycut);
	}

	if(!bDebye && !bDirect){
		Rho_ex->selectPadding(exPadding);
		Rho_ex->setPadding(MyIn.gfftPadding().getMapResidue());
	}
	MySaxs->Allocate(nnx,nny,nnz);
	if(bnoSplineOut)
		MySaxs->SetSplineout();
	Rcut=MyIn.gRcut();
	Rcut_in=MyIn.gRcut_in();
	Clustering=MyIn.bbClust();
	PercolationD::setPercoCutoff(MyIn.gPercoCutoff());

#ifdef HAVE_OPENMP
	const int nthread{omp_get_num_threads()};
	fftw_init_threads();
	fftw_plan_with_nthreads(nthread);
#endif
}

void ExecbSaxsTraj::__SuperCell(MAtoms * atm){
	Matrix co{MetricD{atm->getMt()}.getCO()};
	double mydouble[]={co[XX][XX],co[YY][YY],co[ZZ][ZZ]};
	MySaxs->setSuperCell0(SuperCell0);
	SuperCell=SuperCell0*(*std::max_element(mydouble,mydouble+3));
	try{
		if(SuperCell < 0) throw string("pdb box does not have CRYST1 keyword set, cannot compute.");
	}catch(const string & s){cout << s <<endl;Finale::Finalize::Final();}
	MySaxs->setSuperCell(SuperCell);
}
void ExecbSaxsTraj::__AllocateRho(MAtoms * atm){
	Matrix co{MetricD{atm->getMt()}.getCO()};
	double scx=fabs(SuperCell-co[XX][XX]) < eps?1.0:co[XX][XX]/SuperCell;
	double scy=fabs(SuperCell-co[YY][YY]) < eps?1.0:co[YY][YY]/SuperCell;
	double scz=fabs(SuperCell-co[ZZ][ZZ]) < eps?1.0:co[ZZ][ZZ]/SuperCell;
	int mx=nnx*scx;
	int my=nny*scy;
	int mz=nnz*scz;

	Rho_ex->Allocate(mx,my,mz);
}
void ExecbSaxsTraj::__Compute(const MAtoms * atm){
	switch(WhatToDo){
	case Enums::SAXS:
		MySaxs->ComputeSAXS(Rho_ex,atm);
		break;
	case Enums::SANS:
		MySaxs->ComputeSANS(Rho_ex,atm);
		break;
	case Enums::SQ:
		MySaxs->ComputeSq(Rho_ex,atm);
		break;
	case Enums::ELDENS:
		MySaxs->ComputeSAXS(Rho_ex,atm,true);
		break;


	}
}
void ExecbSaxsTraj::__RunPDB(MAtoms * atm){
	vector<string> data;
	fpdb->clear();
	fpdb->seekg(0);
	for(string str;getline(*fpdb,str);){
		data.push_back(str);
	}
	atm->pdb(data);


	if(Rho_ex) {
		__SuperCell(atm);
		__AllocateRho(atm);
	}

	MySaxs->Setup(atm->getIndx(),Top->get_atSFactor(),WhatToDo==Enums::SANS,WhatToDo==Enums::ELDENS);
	ContactsD * Con0;
	if(Rcut_in < 0) Rcut_in=15.0;
	Con0=new ContactsD(Rcut,Rcut_in);
	Rcut_in=min3(MetricD{atm->getMt()}.getCO());
	Con0->setR(Rcut_in,Rcut_in);
	atm->setrd(*Top);
	if(Clustering){
		atm->SetupPercolate<Enums::noJSON>(*Top);
		atm->Percolate();
		atm->Reconstruct(Con0);
	}
	__Compute(atm);
	MySaxs->Reduce(CurrMPI);
	MySaxs->Averages();
}
void ExecbSaxsTraj::__RunPDBtest(MAtoms * atm){
	vector<string> data;
	fpdb->clear();
	fpdb->seekg(0);
	for(string str;getline(*fpdb,str);){
		data.push_back(str);
	}
	atm->pdb(data);
}
void ExecbSaxsTraj::__RunTrajectory(MAtoms * atm){
	MySaxs->Setup(atm->getIndx(),Top->get_atSFactor(),WhatToDo==Enums::SANS,WhatToDo==Enums::ELDENS);
	myiterators::IteratorAtoms<double> iter_atm(atm,finx,nstart,nend,nskip);
	int MPIsize=CurrMPI->Get_Size();
	int * i_recv=new int [MPIsize];
	double * d_recv=new double [MPIsize];
	ContactsD * Con0;
	if(Rcut_in < 0) Rcut_in=15.0;
	Con0=new ContactsD(Rcut,Rcut_in);

	while((++iter_atm).isReferenced()){
		MAtoms * atmA=iter_atm();

		Con0->setR(Rcut_in,Rcut_in);
		__SuperCell(atmA);
		__AllocateRho(atmA);

		if(Clustering){
			atmA->setrd(*Top);
			static struct Once{Once(MAtoms * atmA,Topol_NS::Topol * myTop){atmA->SetupPercolate();}} _Once(atmA,Top);
			if(bOnce){
				static struct Once_p{explicit Once_p(MAtoms *atmA){atmA->Percolate();}} __Once_p(atmA);
			}else atmA->Percolate();
			atmA->Reconstruct(Con0);
			atmA->CompCM();
		}

		__Compute(atmA);
		int nStep=iter_atm.getTime();
		double Step=atmA->getTime();
		if(MPIsize >1){
			CurrMPI->AllGather(1,&nStep,i_recv);
			CurrMPI->AllGather(1,&Step,d_recv);
			for(auto o=0;o<MPIsize;o++){
				cout << fixed << setw(5) << "----> Step No. " << setw(5) << right << fixed << i_recv[o]
															 << " -- Time = " << d_recv[o] <<" ps \n";
			}
		}else{
			cout << fixed << setw(5) << "----> Step No. " << setw(5) << right << fixed << nStep
														 << " -- Time = " << Step <<" ps \n";
			}
		Rho_ex->Deallocate();
	}
	MySaxs->Reduce(CurrMPI);
	MySaxs->Averages();

	delete [] i_recv;
	delete [] d_recv;
}
void ExecbSaxsTraj::bPrint(ofstream & fout_bin){
	if(!CurrMPI->Get_Rank()){
		MySaxs->bPrintIq(fout_bin);
		cout << "\nProgram completed: Binary output data written to " + fileout_bin << "\n\n";
	}
	CurrMPI->Barrier();
}

ofstream & operator<<(ofstream & fout, ExecbSaxsTraj & y){
	if(!y.CurrMPI->Get_Rank()){
		fout << *y.MySaxs;
		cout << "\nProgram completed: Output data written to " + y.fileout << "\n\n";
	}
	return fout;
}

