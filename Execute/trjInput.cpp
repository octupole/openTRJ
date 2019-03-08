/*
 * trjInput.cpp
 *
 *  Created on: May 22, 2012
 *      Author: marchi
 */

#include "../Execute/trjInput.h"

namespace trj {
trjInput::trjInput(int ntot,char ** v, ClearUsage & clr) {
	auto ClrU=[this](std::initializer_list<int> list){for(auto op: list) Usage[op].clear();};
	vector<string> in;

	inmap["-o"]=in;
	inmap["-obin"]=in;
	inmap["-opdb"]=in;
	inmap["-in"]=in;
	inmap["-dcd"]=in;
	inmap["-xtc"]=in;
	inmap["-b"]=in;
	inmap["-e"]=in;
	inmap["-skip"]=in;
	inmap["-pdb"]=in;
	inmap["-pvol"]=in;
	inmap["-nohyd"]=in;
	inmap["-hyd"]=in;
	inmap["-nodel"]=in;
	inmap["-del"]=in;
	inmap["-rd"]=in;
	inmap["-nord"]=in;
	inmap["-test"]=in;
	inmap["-help"]=in;
	inmap["-parea"]=in;
	inmap["-pgroup"]=in;
	inmap["-surface"]=in;
	inmap["-protein"]=in;
	inmap["-water"]=in;
	inmap["-shell"]=in;
	inmap["-channels"]=in;
	inmap["-solute"]=in;
	inmap["-select"]=in;
	inmap["-detP"]=in;
	inmap["-detO"]=in;
	inmap["-det"]=in;
	inmap["-clust"]=in;
	inmap["-json"]=in;
	inmap["-avgPDB"]=in;
	inmap["-once"]=in;
	inmap["-rho"]=in;
	inmap["-gyro"]=in;

	map<string,vector<string> >::iterator it=inmap.begin();
	for(int n=0;it!=inmap.end();++it,n++){
		Usage[n]=" ";
	}
	int N{0};
	Usage[N++]="\t -dcd file.dcd \n";
	Usage[N++]="\t -xtc file.xtc\n";
	Usage[N++]="\t -pdb file.pdb \n";
	Usage[N++]="\t -o fileout \n"
			"\tFor a file with a .pdb extension: Write pdb files for the solute \n"
			"\tFor a file with a .ndx extention: Write ndx files for the solute \n";
	Usage[N++]="\t -b No. First Frame \n";
	Usage[N++]="\t -e No. Last Frame\n";
	Usage[N++]="\t -in file.bin\n";

	Usage[N++]="\t -skip No. skipped frames \n";
	Usage[N++]="\t -nohyd // Do not include hydrogens in Voronoi [default] \n";
	Usage[N++]="\t -hyd // Include hydrogens in Voronoi\n";
	Usage[N++]="\t -nodel // Do not delete temporary files in MPI runs \n";
	Usage[N++]="\t -del // Delete temporary files in MPI runs [default]\n";
	Usage[N++]="\t -rd // Use atomic radii in the Voronoi calculation [default]\n";
	Usage[N++]="\t -nord // Do not use atomic radii in the Voronoi calculation \n";
	Usage[N++]="\t -solute <string selection>\n"
			"\t\tDefine the residue of the solute to be used in the Voronoi analysis.\n";

	Usage[N++]="\t -test // compute total volume from Voronoi or simulation cell parameters\n";
	Usage[N++]="\t -help // write some on line help \n";
	Usage[N++]="\t -pvol // Print volume \n";
	Usage[N++]="\t -parea // Print area \n";
	Usage[N++]="\t -det filename// Define the polar segment of a detergent residue \n"
			"\t\tThe file filename contains two strings:\n"
			"\t\t   i) the first string is the name of the detergent residue \n"
			"\t\t   ii)the second string contains the polar atoms of the residue\n"
			"\t\t*Important*: The remaining atoms of the residues are considered as hydrophobic\n";
	Usage[N++]="\t -detP <string>// Define a name for the polar segment of a detergent residue \n";
	Usage[N++]="\t -detO <string>// Define a name for the hydrophobic segment of a detergent residue \n";
	Usage[N++]="\t -shell [2]<int n>// Compute number and average volume of the n-shell water\n"
			"\t\tneighbours around solute \n";
	Usage[N++]="\t -clust <float cutoff [0.0] in \u00C5ngstr\u00F6ms>\n"
			"\t\tDo clustering by percolation of the solute. This is used to compute the surfaces and volumes\n"
			"\t\tof each cluster. If the cutoff is zero, the default, the cutoff is chosen according to the\n"
			"\t\tLennard-Jones sigma parameters multiplied by an offset factor of 1.5\n";
	Usage[N++]="\t -json fileout \n"
			"\t\tWrite output to a JSON format\n";
	Usage[N++]="\t -gyro \n"
			"\t\tCompute gyration radius for the molecular clusters\n";
	Usage[N++]="\t -avgPDB Write averaged pdb \n";
	Usage[N++]="\t -once Do percolation only at the beginning \n";
	Usage[N++]="\t -rho <string filename> <float cutoff> <float dx> Compute density solvent profile around cluster \n";
	clr.ClrU(Usage);
	int n=1;
	string key;
	for(;n<ntot;){
		string tmp0(v[n]);
		if(v[n][0] =='-'){
			key.assign(v[n]);
			if(inmap.find(key) != inmap.end()){
				if(inmap[key].empty()) inmap[key].push_back(tmp0);
			} else{
				unknownC.push_back(key);
				inmap[key].push_back(tmp0);
			}
		}
		else{
			inmap[key].push_back(tmp0);
		}
		n++;
	}

}

vector<string> trjInput::getUsage(){
	vector<string> use;
	map<int,string>::iterator it=Usage.begin();
	for(;it!=Usage.end();++it){
		use.push_back(it->second);
	}
	return use;
}
trjInput::~trjInput() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */
