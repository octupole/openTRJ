/*
 * VoronoiMicelles.cpp
 *
 *  Created on: Nov 19, 2017
 *      Author: marchi
 */

#include <VoronoiMicelles.h>

namespace Voro{
const double SURFFACTOR{0.08};

VoronoiMicelles::VoronoiMicelles(Topol & myTop, bool bH){
	nr=myTop.Size();
	Vol=vector<double>(nr);
	Neighs=vector<vector<int>>(nr);
	Surface=vector<vector<double>>(nr);
	label=vector<string>(nr);
	types=vector<int>(nr);
	atTypes=vector<int>(nr);
	SelectedResidues=myTop.gReferenceResidues();


	label=myTop.getAtomName();
	vector<vector<int> > & Index=myTop.gCIndex();
	CIndex=vector<vector<int>>(Index.size());
	for(size_t o=0;o<Index.size();o++){
		for(size_t p=0;p<Index[o].size();p++){
			size_t i=Index[o][p];
			Vol[i]=0.0;
			if(!bH){
				string sub1=label[i];
				if(sub1.find_first_not_of("1234") == 1){
					if(sub1.compare(1,1,"H") != 0) {
						cindex.push_back(i);
						CIndex[o].push_back(i);
					}
				} else{
					if(sub1.compare(0,1,"H") != 0) {
						cindex.push_back(i);
						CIndex[o].push_back(i);
					}
				}
			} else {
				cindex.push_back(i);
				CIndex[o].push_back(i);
			}
		}
	}
	if(bH) {
		cout << std::left << std::fixed<< "\n              " <<std::setw(10) <<
				" Voronoi's calculation includes hydrogens " << endl;
	} else{
		cout << std::left << std::fixed<< "\n              " <<std::setw(10) <<
				" Voronoi's calculation does not includes hydrogens " << endl;
	}
	cout <<"\n" << std::left << std::fixed<<"                  "<<std::setw(10) << " Atoms are  = " << "      "
			<< std::setw(12) << cindex.size() << "\n" <<endl;

	atTypes=myTop.gatResType();
	nresid=myTop.ResSize();

	TypesName=myTop.getTypeNames();
	const vector<int> & tmp=myTop.getAtomTypeNo();
	for(int o=0;o<nr;o++)
		types[o]=tmp[o];
	Residue=myTop.getResinfo();
	Rdii=myTop.getrd();
	nc=ResidueTypes::Size();
	area.Allocate(nresid,nc);
	Vols.Allocate(nresid);
	set<int> tmp1;
	for(string it: this->TypesName){
		tmp1.insert(ResidueTypes::find(it));
	}
	for(auto it:tmp1){
		this->typesResidueMask.push_back(it);
	}

	wShells=vector<vector<int>>(VoronoiSetter::maxLevel);
	vector<int> allAtoms;
	for(auto it0: CIndex)
		for(auto it1: it0)
			allAtoms.push_back(it1);
	this->nWaters=std::count_if(allAtoms.begin(),allAtoms.end(),[this](int i){return types[i] == Water;});
}

void VoronoiMicelles::getData(){

	area=0.0;
	for(int o=0;o<nresid;o++) {
		vector<int> & cindex=CIndex[o];
		double sum_v=0.0;
		vector<tArea> cdx(cindex.begin(),cindex.end());
		std::sort(cdx.begin(),cdx.end(),tAcomp());
		for(unsigned int ia=0;ia<cdx.size();ia++){
			int i=cdx[ia].n;
			vector<tArea> nei0,nei;
			sum_v+=Vol[i];
			for(unsigned int p=0;p<Neighs[i].size();p++)
				nei0.push_back(tArea(Neighs[i][p],Surface[i][p]));

			sort(nei0.begin(),nei0.end(),tAcomp());
			set_difference(nei0.begin(),nei0.end(),cdx.begin(),cdx.end(),back_inserter(nei),tAcomp());

			vector<tArea>::iterator it=nei.begin();
			for(;it != nei.end(); ++it)
				area[o][types[it->n]]+=it->a;

		}
		Vols[o]=sum_v;
	}
	NeighsRd=vector<vector<int>> (Neighs.size());
	for(size_t o{0};o< Neighs.size();o++){
		auto totS=std::accumulate(Surface[o].begin(),Surface[o].end(),0.0);
		for(size_t p{0};p< Surface[o].size();p++)
			if(Surface[o][p] > totS*SURFFACTOR) NeighsRd[o].push_back(Neighs[o][p]);
	}
	Neighsx=&NeighsRd;
	if(VoronoiSetter::bPrintShell) __compShell();
	if(!VolClusters.empty()) __computeAggregate();
}


void VoronoiMicelles::__computeAggregate(){
	for(size_t o=0;o<this->Clusters.size() ;o++) {
		vector<int> & cindex=this->Clusters[o];
		for(unsigned int ia=0;ia<cindex.size();ia++){
			int n=cindex[ia];
			int atCluster_o=atClusters[n];
			VolClusters[o]+=Vol[n];

			for(unsigned int p=0;p<Neighs[n].size();p++){
				int m=Neighs[n][p];
				if(atClusters[m] == atCluster_o) continue;
				int type_n=types[m];
				SurfaceClusters[o][type_n]+=Surface[n][p];
				AreaClusters[o]+=Surface[n][p];
			}

		}
	}
}
void VoronoiMicelles::__searchNeighs(int Level,int n){
	Level++;
	if(Level > VoronoiSetter::maxLevel) return;
	for(unsigned int p=0;p<(*Neighsx)[n].size();p++){
		int o=(*Neighsx)[n][p];
		if(types[o] == Water) {
			wShells[Level-1].push_back(o);
			__searchNeighs(Level,o);
		}
	}
};

void VoronoiMicelles::__compShell(){
	for(auto & it: wShells) it.clear();
	for(size_t o=0;o<CIndex.size() ;o++) {
		vector<int> & cindex=CIndex[o];
		for(unsigned int ia=0;ia<cindex.size();ia++){
			int n=cindex[ia];
			int type_n=types[n];
			if(type_n == Water) continue;
			int Level=0;
			__searchNeighs(Level,n);
		}
	}
	int mWaters{0};
	for(int o{0};o< VoronoiSetter::maxLevel;o++){
		std::sort(wShells[o].begin(),wShells[o].end());
		auto it = std::unique (wShells[o].begin(), wShells[o].end());
		wShells[o].resize( std::distance(wShells[o].begin(),it) );
		for(int p{o-1};p>= 0;p--){
			auto it=std::set_difference(wShells[o].begin(),wShells[o].end(),wShells[p].begin()
					,wShells[p].end(),wShells[o].begin());
			wShells[o].resize(it-wShells[o].begin());
		}

		if(mWaters+=wShells[o].size() == nWaters) break;
	}
}


void VoronoiMicelles::WriteIt(std::ofstream & fout){
	static bool firsttime{true};
	if(firsttime){
		fout << "#  Defined Residue Types: "<<endl;
		fout << "#                 ";
		for(int o=0;o<nc;o++) {
			if(nc-o == 1)
				fout  <<ResidueTypes::getType(o)<< "("<< setw(1)<< o <<") ";
			else
				fout  <<ResidueTypes::getType(o)<< "("<< setw(1)<< o <<"), ";
		}
		fout << "#                 ";
		fout << endl;
		firsttime=false;
	}
	fout << endl;
	fout << "######>> At step No. " << setw(10) << setprecision(2) << fixed<< time << endl;
	int po=0;
	for(size_t o0=0;o0<this->SelectedResidues.size() ;o0++) {
		int o=this->SelectedResidues[o0];
		if(!VoronoiSetter::bPrintVols)continue;
		if(getTypesRes(o) != VoronoiSetter::pGroup && VoronoiSetter::pGroup != -1) continue;
		double a=Vols[o]*1000.0;
		string l=Residue[o];
		(po%4)?fout << setw(10) << setprecision(4) << fixed<<  a << ' ' << setw(4) << l << ' ' << setw(5) << o+1<< ' ':
		  fout << endl << "%$VolRes " << setw(10) << setprecision(4) << fixed << a << ' ' << setw(4) << l  <<' ' << setw(5) << o+1 << ' ';
		po++;
	}
	fout << endl;
	if(!VolClusters.empty()){
		int po{0};
		fout <<endl;
		fout << "#  Volumes of the cluster aggregates: "<<endl;
		fout << "#  No of Clusters = " << fixed << left << VolClusters.size()<< " " <<endl;
		fout << "#  " ;
		for(size_t o{0};o< VolClusters.size();o++){
			if(po%8){
				fout << setw(2)  << fixed<<  o+1 << ' ';
				fout << setw(10) << fixed<<setprecision(2) << VolClusters[o]*1000.0 << " ";
			} else{
			  fout << endl << "%$VolClust " << setw(2) << fixed << o+1 << ' ' ;
			  fout << setw(10) << fixed << setprecision(2) << VolClusters[o]*1000.0 << ' ';
			}
			po++;

		}
		fout << endl;
		fout << endl;

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
		fout<<endl;
		fout << setw(10) << fixed<<left << "#wShell ";
		fout << setw(10) << fixed << "ShellNo";
		fout << setw(10) << fixed << "NoWat" ;
		fout << setw(10) << setprecision(2) << fixed << "Vol/water" << " \n";
		for(size_t o{0};o< wShells.size();o++){
			double vol=Vol0[o]*1000.0;
			fout << setw(10) <<left << fixed << " ";
			fout << setw(10) <<left << fixed <<o;
			fout << setw(10) <<left << fixed << wShells[o].size() ;
			fout << setw(10) <<left << setprecision(2) << fixed
					<<(wShells[o].size()!=0?vol/static_cast<double>(wShells[o].size()):0.0) << " \n";
		}

		fout<<endl;
	}
	if(VoronoiSetter::bPrintAreas ){
		fout << "# Area format:            ";
		for(size_t o0=0;o0<this->typesResidueMask.size() ;o0++) {
			int o=this->typesResidueMask[o0];
			stringstream ss;
			ss<<o;
			fout << setw(1) << fixed << right<<setw(12)<<ResidueTypes::getType(o)+"(" + ss.str()+ ") ";
		}

		fout << endl;

		for(size_t o0=0;o0<SelectedResidues.size() ;o0++){
			fout  << "%$AreaRes " ;
			int o{SelectedResidues[o0]};
			string l=Residue[o];
			int o_type=ResidueTypes::find(l);
			if(VoronoiSetter::pGroup != -1 && o_type != VoronoiSetter::pGroup) continue;
			fout  << right << setw(5)<< l << " " << setw(5) << fixed << ResidueTypes::getType(o_type) << ' ' << setw(3) << fixed << o+1 << ' ' ;
			for(size_t p0=0;p0<this->typesResidueMask.size();p0++) {
				int p=this->typesResidueMask[p0];
				double a=area[o][p]*100.0;
				fout << setw(11) << setprecision(5) << fixed << a << ' ';
			}
			fout << endl;
		}
		fout<<endl;
		if(!AreaClusters.empty()){
			fout<<endl;

			fout << "#                           "<<endl;
			fout << "# Cluster Area format:      ";
			for(size_t o0=0;o0<this->typesResidueMask.size() ;o0++) {
				int o=this->typesResidueMask[o0];
				stringstream ss;
				ss<<o;
				fout << fixed <<setw(11)<<ResidueTypes::getType(o)+"("+ss.str()+") ";
			}
			fout << "     Total "<<endl;
			fout << "# "<<endl;
			for(size_t o{0};o<SurfaceClusters.size();o++){
				fout  << "%$AreaClust              " ;
				for(size_t p0=0;p0<this->typesResidueMask.size();p0++) {
					int p=this->typesResidueMask[p0];
					double a=SurfaceClusters[o][p]*100.0;
					fout << setw(11) << setprecision(2) << right<<fixed << a << ' ';
				}
				double a=AreaClusters[o]*100.0;
				fout << setw(11) << setprecision(2) << fixed <<right << a << ' ';
				fout << endl;

			}
			fout<<"\n\n";
		}
	}

	fout << "# Volume of selection " << endl <<"#        ";
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		fout << fixed<<setw(12)<<right<<ResidueTypes::getType(o)<<"(" << setw(1)<< o << ")";
	}

	fout << "\n"<< endl;
	fout << "%$TotVol ";
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		fout << setw(15) << setprecision(4) << fixed << VolSel[o];
	}
	fout << endl;

	fout << "\n"<< endl;
	fout << "# Interface area of selection " << endl << "#           ";
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		fout << fixed<<setw(10)<<right<<ResidueTypes::getType(o)<<"(" << setw(1)<< o << ")";
	}
	fout << "\n"<< endl;
	for(size_t o0{0};o0<this->typesResidueMask.size();o0++) {
		int o=this->typesResidueMask[o0];
		fout << "#   ";
		fout <<setw(5) <<fixed << left<<ResidueTypes::getType(o)<<"(" << setw(1)<< o << ")";
		for(size_t p0{0};p0<o0;p0++) {
			fout <<setw(13)  << " ";
		}
		for(size_t p0{o0};p0<this->typesResidueMask.size();p0++) {
			int p=this->typesResidueMask[p0];
			fout << right<<setw(13) << interface[o][p] ;
		}
		fout << endl;
	}
	fout << "\n"<< endl;

}

VoronoiMicelles::~VoronoiMicelles() {
	// TODO Auto-generated destructor stub
}

}
