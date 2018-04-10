/*
 * TrjRead.cpp
 *
 *  Created on: Dec 23, 2015
 *      Author: marchi
 */

#include "TrjReadSaxs.h"
namespace trj {
Parallel::NewMPI * TrjRead::CurrMPI=nullptr;

TrjRead::TrjRead(int nv,char ** v): trjInput::trjInput(nv,v) {
	// TODO Auto-generated constructor stub
	string comm0(v[0]);
	size_t mypos=comm0.find_last_of("/")+1;
	size_t length=comm0.size()-mypos;
	string command=comm0.substr(mypos,length);
	string errmsg;
	string Usage="Usage:\t"+ command + "\n";
	vector<string> use=this->getUsage();
	vector<string> SelRes;
	string Reference;
	for(unsigned int n=0;n<use.size();n++)
		Usage+=use[n];
	Usage+="\n\t Default values in square brackets []\n";
	try{
		if(nv == 1) throw Usage;
		if(int m=this->bTest().size()) {
			errmsg=" Command(s) not found: ";
			for(unsigned int n=0;n<m;n++)
				errmsg+=this->bTest()[n]+"  ";
			errmsg+="\n"+Usage;
			throw errmsg;
			}
		}
	catch(const string & s){
		cout << s << endl;
		Finale::Finalize::Final();
	}

}
void TrjRead::Input(){
	ifstream ftest,fdefdomain;
#ifdef __INTEL  //  set the conversion factors by hand as until intel c++ compiler version 16 user-defined literals are not allowed
	struct dummy{
		double val;
	};
	dummy e0_conversion{0.694461547828659};
	ConvFactor=55.7036562713774;
#else
	auto bimbo=1.0_el/(1.0_nm*1.0_nm);

	auto fluct_en=(1.0/(4.0*M_PI))*bimbo*1.0_el*1.0_nm/eps0;
	auto Factor=fluct_en/(kt300/avogad) ;
	ConvFactor=Factor.val;

	auto FieldConv=1.0_Volt/1.0_nm;
	auto e0_conversion=4.0*M_PI*FieldConv*eps0/bimbo;
#endif
	try{
		if(!inmap["-dcd"].empty()) {
			if(inmap["-dcd"].size() < 2) throw string("\n Filename expected for " + inmap["-dcd"][0] + " option \n");
			if(inmap["-dcd"].size() > 2) throw string("\n More than one entry for " + inmap["-dcd"][0] + " option \n");
			filein=inmap["-dcd"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
			ftest.close();

			finx=new FstreamF(filein);
			inputfile=true;
		}
		if(!inmap["-in"].empty()) {
			if(inmap["-in"].size() < 2) throw string("\n at least one filename expected for " + inmap["-in"][0] + " option\n ");
			if(inmap["-in"].size() > 3) throw string("\n More than two entries for " + inmap["-in"][0] + " option \n");
			filein=inmap["-in"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open system " + filein + "!!\n");
			ftest.close();
			fin1=new ifstream(filein,ios_base::binary);
			inputfile=true;
			bPost=true;
			if(inmap["-in"].size() == 3){
				filein=inmap["-in"][2];
				ftest.open(filein.c_str(),ios::in);
				if(!ftest) throw string("\n Cannot open contrast " + filein + "!!\n");
				ftest.close();
				fin2=new ifstream(filein,ios_base::binary);
			}
		}
		if(!inmap["-padd"].empty()) {
			if(inmap["-padd"].size() < 2) throw string("\n filename expected for " + inmap["-padd"][0] + " option\n ");
			if(inmap["-padd"].size() > 2) throw string("\n More than one entry for " + inmap["-padd"][0] + " option \n");
			string myParam=inmap["-padd"][1];
			if(myParam == "pbc") {
				WhichPadding=Enums::Periodic;
			} else if(myParam == "zero")
				WhichPadding=Enums::zero;
			else if(myParam == "avg")
				WhichPadding=Enums::avgDensity;
			else{
				WhichPadding=Enums::myDensity;
				filein=inmap["-padd"][1];
				ftest.open(filein.c_str(),ios::in);
				if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
				ftest.close();
				fin_padding=new ifstream(filein,ios_base::binary);
			}
		}
		if(!inmap["-xtc"].empty()) {
			if(inmap["-xtc"].size() < 2) throw string("\n filename expected for " + inmap["-xtc"][0] + " option\n ");
			if(inmap["-xtc"].size() > 2) throw string("\n More than one entry for " + inmap["-xtc"][0] + " option \n");
			filein=inmap["-xtc"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
			ftest.close();
			finx=new FstreamC(filein,"rb");
			inputfile=true;
		}
		if(!inmap["-xtc2"].empty()) {
			if(inmap["-xtc2"].size() < 2) throw string("\n filename expected for " + inmap["-xtc2"][0] + " option\n ");
			if(inmap["-xtc2"].size() > 2) throw string("\n More than one entry for " + inmap["-xtc2"][0] + " option \n");
			filein=inmap["-xtc2"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
			ftest.close();
			fin2x=new FstreamC(filein,"rb");
			inputfile=true;
		}
		if(!inmap["-rcm"].empty()) {
			if(inmap["-rcm"].size() < 2) throw string("\n filename expected for " + inmap["-rcm"][0] + " option\n ");
			if(inmap["-rcm"].size() > 2) throw string("\n More than one entry for " + inmap["-rcm"][0] + " option \n");
			filein=inmap["-rcm"][1];
			ftest.open(filein.c_str(),ios::in);
			if(!ftest) throw string("\n Cannot open " + filein + "!!\n");
			ftest.close();
			fin_rcmx=new ifstream(filein,ios::in | ios::binary);
			inputfile=true;
		}
		if(!inmap["-pdb"].empty()) {
			if(inmap["-pdb"].size() < 2) throw string("\n filename expected for " + inmap["-pdb"][0] + " option \n");
			if(inmap["-pdb"].size() > 2) throw string("\n More than one entry for " + inmap["-pdb"][0] + " option \n");
			filepdb=inmap["-pdb"][1];
			fpdb.open(filepdb.c_str(),ios::in);
			if(!fpdb) throw string("\n Cannot open " + filepdb + "!!\n");
		}
		if(!inmap["-select"].empty()) {
			if(inmap["-select"].size() != 2) throw string(" String of selected residues needed for " + inmap["-select"][0] + " option ");
			string selection=inmap["-select"][1];
			stringstream iss(selection);
			copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter<vector<string> >(SelRes));
		}
		if(!inmap["-solute"].empty()) {
			if(inmap["-solute"].size() != 2) throw string(" Reference residue for Micelle centering needed " + inmap["-solute"][0] + " option ");
			string Ref=inmap["-solute"][1];
			stringstream iss(Ref);
			copy(istream_iterator<string>(iss),
					istream_iterator<string>(),
					back_inserter<vector<string> >(Reference));
		}
		if(!inmap["-o"].empty()) {
			if(inmap["-o"].size() < 2) throw string("\n filename expected for " + inmap["-o"][0] + " option \n");
			if(inmap["-o"].size() > 2) throw string("\n More than one entry for " + inmap["-o"][0] + " option \n");
			fileout=inmap["-o"][1];
		}
		if(!inmap["-obin"].empty()) {
			if(inmap["-obin"].size() < 2) throw string("\n filename expected for " + inmap["-obin"][0] + " option \n");
			if(inmap["-obin"].size() > 2) throw string("\n More than one entry for " + inmap["-obin"][0] + " option \n");
			fileout_bin=inmap["-obin"][1];
			bOutBin=true;
		}
		if(!inmap["-ic"].empty()) {
			if(inmap["-ic"].size() < 2) throw string("\n filename expected for " + inmap["-ic"][0] + " option \n");
			if(inmap["-ic"].size() > 2) throw string("\n More than one entry for " + inmap["-ic"][0] + " option \n");
			string fileContrast=inmap["-ic"][1];
			bContrast=true;
		}
		if(!inmap["-chg"].empty()) {
			if(inmap["-chg"].size() < 2) throw string("\n filename expected for " + inmap["-chg"][0] + " option \n");
			if(inmap["-chg"].size() > 2) throw string("\n More than one entry for " + inmap["-chg"][0] + " option \n");
			filechg=inmap["-chg"][1];
			bdip=true;
		}
		if(!inmap["-dip"].empty()) {
			if(inmap["-dip"].size() != 1) throw string(" No argument for " + inmap["-dip"][0] + " option ");
			bdip=true;
		}
		if(!inmap["-b"].empty()) {
			if(inmap["-b"].size() != 2) throw string(" Number of first frame needed for " + inmap["-b"][0] + " option ");
			stringstream(inmap["-b"][1])>> nstart;
		}
		if(!inmap["-diffk"].empty()) {
			if(inmap["-diffk"].size() != 1) throw string(" No argument for " + inmap["-diffk"][0] + " option ");
			WhichDiffusion=diffk;

		}
		if(!inmap["-diffq"].empty()) {
			if(inmap["-diffq"].size() != 2) throw string(" Three Wigner arguments in a string are needed for " + inmap["-diffq"][0] + " option ");
			WhichDiffusion=diffq;
			stringstream(inmap["-diffq"][1])>> WignerArgs[0]>> WignerArgs[1] >> WignerArgs[2];
		}
		if(!inmap["-diffm"].empty()) {
			if(inmap["-diffm"].size() != 2) throw string(" Three Wigner arguments in a string are needed for " + inmap["-diffm"][0] + " option ");
			WhichDiffusion=diffMicelles;
			stringstream(inmap["-diffm"][1])>> WignerArgs[0]>> WignerArgs[1] >> WignerArgs[2];
		}
		if(!inmap["-e"].empty()) {
			if(inmap["-e"].size() != 2) throw string(" Number of last frame needed for " + inmap["-e"][0] + " option ");
			stringstream(inmap["-e"][1])>> nend;
		}
		if(!inmap["-skip"].empty()) {
			if(inmap["-skip"].size() != 2) throw string(" Number of skipped frames needed for " + inmap["-skip"][0] + " option ");
			stringstream(inmap["-skip"][1])>> nskip;
		}
		if(!inmap["-nohyd"].empty()) {
			if(inmap["-nohyd"].size() != 1) throw string(" No argument for " + inmap["-nohyd"][0] + " option ");
				bHyd=false;
		}
		if(!inmap["-hyd"].empty()) {
			if(inmap["-hyd"].size() != 1) throw string(" No argument for " + inmap["-hyd"][0] + " option ");
				bHyd=true;
		}
		if(!inmap["-nodel"].empty()) {
			if(inmap["-nodel"].size() != 1) throw string(" No argument for " + inmap["-nodel"][0] + " option ");
				bDel=false;
		}
		if(!inmap["-del"].empty()) {
			if(inmap["-del"].size() != 1) throw string(" No argument for " + inmap["-del"][0] + " option ");
				bDel=false;
		}
		if(!inmap["-neq"].empty()) {
			if(inmap["-neq"].size() != 1) throw string(" No argument for " + inmap["-neq"][0] + " option ");
				bNonEq=true;
		}
		if(!inmap["-prof"].empty()) {
			if(inmap["-prof"].size() != 1) throw string(" No argument for  " + inmap["-prof"][0] + " option ");
			bprof=true;
		}

		if(!inmap["-cut"].empty()) {
			if(inmap["-cut"].size() != 2) throw string(" General cutoff needed for  " + inmap["-cut"][0] + " option ");
			stringstream(inmap["-cut"][1])>> cutYZ;
			cutYZ*=unit_nm;
		}
		if(!inmap["-bsp"].empty()) {
			if(inmap["-bsp"].size() != 1) throw string(" no argument for " + inmap["-bsp"][0] + " option ");
			bBSP=true;
		}
		if(!inmap["-rcut"].empty()) {
			if(inmap["-rcut"].size() != 2) throw string(" Linked list cutoff needed for  " + inmap["-rcut"][0] + " option ");
			stringstream(inmap["-rcut"][1])>> Rcut;
			Rcut*=unit_nm;
		}
		if(!inmap["-molw"].empty()) {
			if(inmap["-molw"].size() != 2) throw string(" Solute molecular weight needed for  " + inmap["-molw"][0] + " option ");
			stringstream(inmap["-molw"][1])>> MassSolute;
		}
		if(!inmap["-rin"].empty()) {
			if(inmap["-rin"].size() != 2) throw string(" Linked list cutoff needed for  " + inmap["-rin"][0] + " option ");
			stringstream(inmap["-rin"][1])>> Rcut_in;
			Rcut_in*=unit_nm;
		}

		if(!inmap["-filter"].empty()) {
			if(inmap["-filter"].size() != 2) throw string(" Sigma Gaussian width in Angstroems needed for  " + inmap["-filter"][0] + " option ");
			stringstream(inmap["-filter"][1])>> ewSigma;
			ewSigma*=unit_nm;
		}

		if(!inmap["-cutX"].empty()) {
			if(inmap["-cutX"].size() != 2) throw string(" Out of plane cutoff needed for  " + inmap["-cutX"][0] + " option ");
			stringstream(inmap["-cutX"][1])>> cutX;
			cutX*=unit_nm;
		}
		if(!inmap["-qhist"].empty()) {
			if(inmap["-qhist"].size() != 3) throw string(" Bin size and cutoff in q-space needed for  " + inmap["-qhist"][0] + " option ");
			stringstream(inmap["-qhist"][1])>> Myd;
			stringstream(inmap["-qhist"][2])>> Mycut;
		}
		if(!inmap["-pdx1"].empty()) {
			if(inmap["-pdx1"].size() != 2) throw string(" Bin size needed for  " + inmap["-pdx1"][0] + " option ");
			stringstream(inmap["-pdx1"][1])>> ddx1;
			ddx1*=unit_nm;
		}
		if(!inmap["-pdx2"].empty()) {
			if(inmap["-pdx2"].size() != 2) throw string(" Bin size needed for  " + inmap["-pdx2"][0] + " option ");
			stringstream(inmap["-pdx2"][1])>> ddx2;
			ddx2*=unit_nm;
		}
		if(!inmap["-e0"].empty()) {
			int Nm=inmap["-e0"].size();
			if(Nm == 4){
				stringstream(inmap["-e0"][1])>> e0[XX];
				stringstream(inmap["-e0"][2])>> e0[YY];
				stringstream(inmap["-e0"][3])>> e0[ZZ];
				e0=e0*e0_conversion.val;

			} else throw string(" Vector components needed for  " + inmap["-e0"][0] + " option ");
		}
		if(!inmap["-grid"].empty()) {
			int Nm=inmap["-grid"].size();
			try{
				if(Nm == 2){
					stringstream(inmap["-grid"][1])>> nnx;
					nny=nnx; nnz=nnx;
				} else if(Nm == 4){
					stringstream(inmap["-grid"][1])>> nnx;
					stringstream(inmap["-grid"][2])>> nny;
					stringstream(inmap["-grid"][3])>> nnz;
				} else throw string("\n Warning: No Grid points provided " + inmap["-grid"][0] + " option. "
						+"Hope this is correct \n");
			}catch(const string & s){
				cout << s <<endl;
			}
		}
		if(!inmap["-irl"].empty()) {
			if(inmap["-irl"].size() != 2) throw string(" filename expected for " + inmap["-irl"][0] + " option ");
			string fileirl=inmap["-irl"][1];
			firl=new ifstream;
			firl->open(fileirl.c_str(),ios::in);
			if(!(*firl)) throw string("\n Cannot open " + fileirl + "!!\n");
		}
		if(!inmap["-idb"].empty()) {
			if(inmap["-idb"].size() != 2) throw string(" filename expected for " + inmap["-idn"][0] + " option ");
			string fileidb=inmap["-idb"][1];
			fidb=new ifstream;
			fidb->open(fileidb.c_str(),ios::in);
			if(!(*fidb)) throw string("\n Cannot open " + fileidb + "!!\n");
		}
		if(!inmap["-Voro"].empty()) {
			if(inmap["-Voro"].size() != 1) throw string(" No argument for " + inmap["-Voro"][0] + " option ");
				bVoro=true;
		}
		if(!inmap["-clust"].empty()) {
			if(inmap["-clust"].size() > 2) throw string("\n Only one optional argument is allowed for " + inmap["-clust"][0] + "\n");
			bClust=true;
			if(inmap["-clust"].size() == 2) {
				stringstream(inmap["-clust"][1])>> PercoCutoff;
				PercoCutoff/=10.0;
			}
		}
		if(!inmap["-rho"].empty()) {
			int Nm=inmap["-rho"].size();
			bRho=true;
			switch(Nm){
			case 1:
				break;
			case 2:
				stringstream(inmap["-rho"][1])>> MyOrder;
				break;
			case 3:
				stringstream(inmap["-rho"][1])>> MyOrder;
				stringstream(inmap["-rho"][2])>> ewSigma;
				ewSigma*=unit_nm;
				break;
			default:
				throw string("\n At most 2 arguments for " + inmap["-rho"][0] + " option \n");
			}

		}
		if(!inmap["-lagrange"].empty()) {
			int Nm=inmap["-lagrange"].size();
			bSaxs=true;
			switch(Nm){
			case 1:
				break;
			case 2:
				stringstream(inmap["-lagrange"][1])>> MyOrder;
				break;
			case 3:
				stringstream(inmap["-lagrange"][1])>> MyOrder;
				stringstream(inmap["-lagrange"][2])>> Mycut;
				break;
			case 4:
				stringstream(inmap["-lagrange"][1])>> MyOrder;
				stringstream(inmap["-lagrange"][2])>> Myd;
				stringstream(inmap["-lagrange"][3])>> Mycut;
				break;
			default:
				stringstream ss;
				ss<< Nm;
				throw string("\n At most 3 arguments for " + inmap["-lagrange"][4] + " option, but found "+ ss.str()+ "\n");
			}

		}
		if(!inmap["-spline"].empty()) {
			int Nm=inmap["-spline"].size();
			bSaxs=true;
			bSaxsBSP=true;
			switch(Nm){
			case 1:
				break;
			case 2:
				stringstream(inmap["-spline"][1])>> MyOrder;
				break;
			case 3:
				stringstream(inmap["-spline"][1])>> MyOrder;
				stringstream(inmap["-spline"][2])>> Mycut;
				break;
			case 4:
				stringstream(inmap["-spline"][1])>> MyOrder;
				stringstream(inmap["-spline"][2])>> Myd;
				stringstream(inmap["-spline"][3])>> Mycut;
				break;
			default:
				throw string("\n At most 3 arguments for " + inmap["-spline"][0] + " option \n");
			}

		}

		if(!inmap["-saxs"].empty()) {
			int Nm=inmap["-saxs"].size();
			bSaxs=true;
			switch(Nm){
			case 1:
				break;
			case 2:
				stringstream(inmap["-saxs"][1])>> MyOrder;
				break;
			case 3:
				stringstream(inmap["-saxs"][1])>> MyOrder;
				stringstream(inmap["-saxs"][2])>> Mycut;
				break;
			case 4:
				stringstream(inmap["-saxs"][1])>> MyOrder;
				stringstream(inmap["-saxs"][2])>> Myd;
				stringstream(inmap["-saxs"][3])>> Mycut;
				break;
			default:
				stringstream ss;
				ss<< Nm;
				throw string("\n At most 3 arguments for " + inmap["-saxs"][4] + " option, but found "+ ss.str()+ "\n");
			}

		}
		if(!inmap["-saxssp"].empty()) {
			int Nm=inmap["-saxssp"].size();
			bSaxs=true;
			bSaxsBSP=true;
			switch(Nm){
			case 1:
				break;
			case 2:
				stringstream(inmap["-saxssp"][1])>> MyOrder;
				break;
			case 3:
				stringstream(inmap["-saxssp"][1])>> MyOrder;
				stringstream(inmap["-saxssp"][2])>> Mycut;
				break;
			case 4:
				stringstream(inmap["-saxssp"][1])>> MyOrder;
				stringstream(inmap["-saxssp"][2])>> Myd;
				stringstream(inmap["-saxssp"][3])>> Mycut;
				break;
			default:
				throw string("\n At most 3 arguments for " + inmap["-saxssp"][0] + " option \n");
			}

		}
		if(!inmap["-what"].empty()) {
			bSaxs=true;
			string myWhat=inmap["-what"][1];
			if(myWhat == "sq"){
				bClust=true;
				WhatToDo=Enums::SQ;
			}
			else if(myWhat == "saxs")
				WhatToDo=Enums::SAXS;
			else if(myWhat == "sans")
				WhatToDo=Enums::SANS;
			else{
				throw string("\n -what can only choose argument SQ,SAXS or SANS not "+myWhat+" !! ");
			}
		}
		if(!inmap["-debye"].empty()) {
			int Nm=inmap["-debye"].size();
			bSaxs=true;
			bDebye=true;
			switch(Nm){
			case 1:
				break;
			case 2:
				stringstream(inmap["-debye"][1])>> Mycut;
				break;
			case 3:
				stringstream(inmap["-debye"][1])>> Myd;
				stringstream(inmap["-debye"][2])>> Mycut;
				break;
			default:
				stringstream ss;
				ss<< Nm;
				throw string("\n At most 2 arguments for " + inmap["-debye"][4] + " option, but found "+ ss.str()+ "\n");
			}

		}
		if(!inmap["-direct"].empty()) {
			int Nm=inmap["-direct"].size();
			bSaxs=true;
			bDirect=true;
			switch(Nm){
			case 1:
				break;
			case 2:
				stringstream(inmap["-direct"][1])>> Mycut;
				break;
			case 3:
				stringstream(inmap["-direct"][1])>> Myd;
				stringstream(inmap["-direct"][2])>> Mycut;
				break;
			default:
				stringstream ss;
				ss<< Nm;
				throw string("\n At most 2 arguments for " + inmap["-direct"][4] + " option, but found "+ ss.str()+ "\n");
			}

		}

		if(!inmap["-lowq"].empty()) {
			int Nm=inmap["-lowq"].size();
			if(Nm == 1)
				BoxMultiply=3;
			else if(Nm == 2){
				stringstream(inmap["-lowq"][1])>> BoxMultiply;
			}
			else throw string(" Only one parameters is needed for  " + inmap["-lowq"][0] + " option ");
		}
		if(!inmap["-supcell"].empty()) {
			int Nm=inmap["-supcell"].size();
			if(Nm == 2){
				stringstream(inmap["-supcell"][1])>> SuperCellSide;
			}
			else throw string(" One parameters is needed for  " + inmap["-supcell"][0] + " option ");
		}
		if(!inmap["-nosplout"].empty()) {
			if(inmap["-nosplout"].size() != 1) throw string("\n No argument for " + inmap["-nosplout"][0] + " option \n");
			bnoSplineOut=true;
		}
		if(!inmap["-once"].empty()) {
			if(inmap["-once"].size() != 1) throw string("\n No argument for " + inmap["-once"][0] + " option \n");
			bOnce=true;
		}
		if(!inmap["-nofluct"].empty()) {
			if(inmap["-nofluct"].size() != 2) throw string(" Number of first frame needed for " + inmap["-nofluct"][0] + " option ");
			stringstream(inmap["-nofluct"][1])>> nacc;
			bFluct=false;
		}
		if(!inmap["-nosolv"].empty()) {
			if(inmap["-nosolv"].size() != 1) throw string(" No argument for " + inmap["-nosolv"][0] + " option ");
			bFluct=false;
			bFixed=true;
		}

		if(!inmap["-pVct"].empty()) {
			if(inmap["-pVct"].size() != 2) throw string(" Perpendicular axes needed for  " + inmap["-pVct"][0] + " option ");
			sVect=stringstream(inmap["-pVct"][1]).str();
		}
//		if(!inmap["-help"].empty()) {
//			if(inmap["-help"].size() != 1) throw string(" No argument for " + inmap["-help"][0] + " option ");
//				throw Help;
//		}
	}catch(const string & s){
		cout << s << endl;
		Finale::Finalize::Final();
	}
	if(fileout.substr(fileout.find_first_of(".")+1).find("elc") != string::npos){
		fout_elx=new ofstream(fileout.c_str(),ios_base::binary);
		bOutBin=true;
	} else if(fileout.substr(fileout.find_first_of(".")+1).find("prf") != string::npos){
		bprof=true;
	}  else if(fileout.substr(fileout.find_first_of(".")+1).find("clt") != string::npos){
		fclust=new ofstream(fileout.c_str(),ios::out);
	} else if(fileout.substr(fileout.find_first_of(".")+1).find("rcm") != string::npos){
		fout_cmx=new ofstream(fileout.c_str(),ofstream::out);
	} else if(bOutBin && bdip ){
		fout_elx=new ofstream(fileout.c_str(),ios_base::binary);
	}

	string tmp1=fileout.substr(0,fileout.find_first_of("."));
	string tmp2=".dat";
	if(bprof){
		fileoutp1=tmp1+"-1D"+tmp2;
		fileoutp2=tmp1+"-2D"+tmp2;
		fileoutp3=tmp1+"-3D"+tmp2;
		foutp1=new ofstream(fileoutp1.c_str(),ios::out);
		foutp2=new ofstream(fileoutp2.c_str(),ios::out);
		foutp3=new ofstream(fileoutp3.c_str(),ios::out);
	}
	if(bSaxs){
		foutx=new ofstream(fileout.c_str(),ios::out);
		if(bOutBin) {
			if(CurrMPI->AmI_Parallel()) {
				if(!CurrMPI->Get_Rank())
					BackupFile{fileout_bin};
				CurrMPI->Barrier();
			} else BackupFile{fileout_bin};
			fout_binx=new ofstream(fileout_bin.c_str(),ios::out|ios_base::binary);

		}
	}


	gFinx=finx;
	gFin2x=fin2x;
	gFout_xtcx=fout_xtcx;
	gFout_pdbx=fout_pdbx;
	gFout_dip=fout_dip;
	gFout_cmx=fout_cmx;
	gFout_elx=fout_elx;
	gFout_binx=fout_binx;
	gFoutx=foutx;
	gFoutp1=foutp1;
	gFoutp2=foutp2;
	gFoutp3=foutp3;
	gFclust=fclust;
	gFdomain=&fdomain;
	gFirl=firl;
	gFidb=fidb;
	gFin_rcmx=fin_rcmx;
	gFin1=fin1;
	gFin2=fin2;
	gFin_contrast=fin_contrast;
	gFin_padding=fin_padding;
	gFoutsaxs=fout_saxsx;

	gFpdb=&fpdb;
	gFchg=&fchg;
	gFdefdomain=&fdefdomain;
	gFtest=&ftest;
	gE0()=e0;

	if(fin1 && fin2) bdip=false;
	try{

		if(this->gFinx() && (this->gFin1() || this->gFin2() )) throw string(" No input file(s) given program stops! ");
		if(!Reference.empty()) {
			vector<string> ss;
			std::sort(SelRes.begin(),SelRes.end());std::sort(Reference.begin(),Reference.end());
			std::set_intersection(SelRes.begin(),SelRes.end(),Reference.begin(),Reference.end(),back_inserter(ss));
			if(ss.empty()) throw string("Residue in -solute and -select must intersect");
		}
		if(finx && nstart > nend) throw string("Initial trajectory step is larger than the set final step. Change and rerun. ");
		if(WhichPadding == Enums::Periodic) throw string("Periodic padding not yet implemented.");
	} catch(const string & s){
		cout << s <<endl;
		Finale::Finalize::Final();
	}
	try{
		if(Reference.empty() &&  bClust) throw string("Cannot find clusters if solute is not provided. Clustering is off.");
	} catch(const string & s){
		cout << "\n\nWarning: " + s +"\n\n"<<endl;
		bClust=false;
	}
	if(WhichPadding == Enums::myDensity){
		// read padding file to padding density
		vector<string> Padd;vector<double> Natms;
		for(string str;getline(*fin_padding,str);){
			string tmp1;double tmp2;
			if(str.empty()) continue;
			stringstream ss(str);
			ss>> tmp1;ss>>tmp2;
			Padd.push_back(tmp1);
			Natms.push_back(tmp2/(unit_nm*unit_nm*unit_nm));
		}
		fftPadding(Padd,Natms);
	}
	Mycut/=unit_nm;
	Myd/=unit_nm;
}
TrjRead::~TrjRead() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */
