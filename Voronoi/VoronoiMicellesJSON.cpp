/*
 * VoronoiMicellesJSON.cpp
 *
 *  Created on: Dec 17, 2017
 *      Author: marchi
 */

#include "VoronoiMicellesJSON.h"

namespace Voro {

void VoronoiMicellesJSON::WriteIt(std::ofstream & fout){
	static bool firsttime{true};
	json myJson;

	if(firsttime){
		fout << "{";
		json myTypes;
		for(int o=0;o<nc;o++) {
			myTypes["Types"][ResidueTypes::getType(o)]=o;
		}

		for(size_t o0=0;o0<SelectedResidues.size() ;o0++){
			int o{SelectedResidues[o0]};
			string l=Residue[o];
			int o_type=ResidueTypes::find(l);

			myTypes["Res"][l]={o_type,ResidueTypes::getType(o_type)};
			myTypes["Res"]["List"].push_back(l);
		}

		auto it=myTypes.begin();
		for(;it!=myTypes.end();it++){
			if(it != myTypes.begin()) fout << ",";
			fout << "\""<<it.key()<<"\": ";
			fout << *it;
		}
		fout <<","<<endl;
		firsttime=false;
	}
	json & myTime=myJson[std::to_string(time)];
	json & myVol=myJson[std::to_string(time)]["Vol"];
	json & myArea=myJson[std::to_string(time)]["Area"];
	json & myClust=myJson[std::to_string(time)]["Clust"];
	json & myShell=myJson[std::to_string(time)]["Shell"];
	json & myAClust=myJson[std::to_string(time)]["AClust"];
	json & myGlobal=myJson[std::to_string(time)]["Global"];
	int po=0;
	if(VoronoiSetter::bPrintVols){
		for(size_t o0=0;o0<this->SelectedResidues.size() ;o0++) {
			int o=this->SelectedResidues[o0];
			double a=Vols[o]*1000.0;
			myVol.push_back(a);
		}
	}
	if(!VolClusters.empty()){
		for(size_t o{0};o< VolClusters.size();o++){
			myClust.push_back(VolClusters[o]*1000.0);
		}
	}
	vector<double> VolSel(nc,0.0);
	vector<vector<double>> interface(nc,VolSel);
	for(int o=0;o<nresid ;o++){
		int o_type=ResidueTypes::find(Residue[o]);
		VolSel[o_type]+=Vols[o]*1000.0;
		for(int p=0;p<nc;p++) {
			double a=area[o][p]*100.0;
			interface[o_type][p]+=a;
		}
	}
	if(VoronoiSetter::bPrintShell){
		vector<double> Vol0(wShells.size());
		for(size_t o{0};o< wShells.size();o++){
			Vol0[o]=0;
			for(size_t p0{0};p0<wShells[o].size();p0++){
				int p=wShells[o][p0];
				Vol0[o]+=Vol[p];
			}
		}

		for(size_t o{0};o< wShells.size();o++){
			size_t shellSize=wShells[o].size();
			double vol=shellSize!=0?Vol0[o]*1000.0/static_cast<double>(shellSize):0.0;
			myShell[std::to_string(o)]={{"No", shellSize}};
			myShell[std::to_string(o)]+={"Vol", vol};
		}
	}

	if(VoronoiSetter::bPrintAreas ){
		for(size_t o0=0;o0<SelectedResidues.size() ;o0++){

			int o{SelectedResidues[o0]};
			json tmp;
			for(size_t p0=0;p0<this->typesResidueMask.size();p0++) {
				int p=this->typesResidueMask[p0];
				string type=ResidueTypes::getType(p);
				double a=area[o][p]*100.0;
				if(p0==0){
					tmp={{type, a}};
				} else{
					tmp+={type,a};
				}
			}
			myArea.push_back(tmp);
		}
	}
	if(!AreaClusters.empty()){
		for(size_t o{0};o<SurfaceClusters.size();o++){
			size_t size_clust=Clusters[o].size();
			json tmp;
			for(size_t p0=0;p0<this->typesResidueMask.size();p0++) {
				int p=this->typesResidueMask[p0];
				string l=ResidueTypes::getType(p);
				double a=SurfaceClusters[o][p]*100.0;
				tmp[l]=a;

			}
			double a=AreaClusters[o]*100.0;
			tmp["Tot"]=a;
			tmp["NoAtm"]=size_clust;
			myAClust.push_back(tmp);
		}
	}
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		string l=ResidueTypes::getType(o);
		myGlobal["Vol"][l]=VolSel[o];
	}

	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		string lab_o=ResidueTypes::getType(o);
		for(size_t p0{o0};p0<this->typesResidueMask.size();p0++) {
			int p=this->typesResidueMask[p0];
			string lab_p=ResidueTypes::getType(p);
			myGlobal["Iface"][lab_o][lab_p]=interface[o][p] ;
			myGlobal["Iface"][lab_p][lab_o]=interface[o][p] ;
		}
	}

// Append to JSON file
	auto it=myJson.begin();
	for(;it!=myJson.end();it++){
		fout << "\""<<it.key()<<"\": ";
		fout << *it;
	}
	fout << ",";
}
void VoronoiMicellesJSON::WriteLastJSON(std::ofstream & fout){
	size_t pos=fout.tellp();
	fout.seekp(pos-1);
	fout<<"}"<<endl;
}

VoronoiMicellesJSON::~VoronoiMicellesJSON() {
	// TODO Auto-generated destructor stub
}

} /* namespace Voro */
