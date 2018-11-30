/*
 * trjInput.cpp
 *
 *  Created on: May 22, 2012
 *      Author: marchi
 */

#include "trjInput.h"

namespace trj {
trjInput::trjInput(int ntot,char ** v) {
	auto ClrU=[this](std::initializer_list<int> list){for(auto op: list) Usage[op].clear();};
	vector<string> in;

	inmap["-o"]=in;
	inmap["-obin"]=in;
	inmap["-dcd"]=in;
	inmap["-xtc"]=in;
	inmap["-rcm"]=in;
	inmap["-b"]=in;
	inmap["-e"]=in;
	inmap["-skip"]=in;
	inmap["-pdb"]=in;
	inmap["-nohyd"]=in;
	inmap["-hyd"]=in;
	inmap["-nodel"]=in;
	inmap["-del"]=in;
	inmap["-solute"]=in;
	inmap["-select"]=in;
	inmap["-grid"]=in;
	inmap["-clust"]=in;
	inmap["-in"]=in;
	inmap["-in2"]=in;
	inmap["-file"]=in;
	inmap["-neq"]=in;
	inmap["-dip"]=in;
	inmap["-lagrange"]=in;
	inmap["-spline"]=in;
	inmap["-saxs"]=in;
	inmap["-saxssp"]=in;
	inmap["-what"]=in;
	inmap["-supcell"]=in;
	inmap["-padd"]=in;
	inmap["-qhist"]=in;
	inmap["-debye"]=in;
	inmap["-direct"]=in;
	inmap["-nosplout"]=in;
	inmap["-molw"]=in;
	inmap["-nofluct"]=in;
	inmap["-nosolv"]=in;
	inmap["-once"]=in;
	inmap["-dens"]=in;
	inmap["-h"]=in;

	map<string,vector<string> >::iterator it=inmap.begin();
	size_t i{0};
	Usage[i++]="\t -dcd <string filename>\n"
			"\t\tInput the trajectory file in NAMD format. This is experimental and \n"
			"\t\thas not yet been thoroughly tested.\n";
	Usage[i++]="\t -xtc <string filename> \n"
			"\t\tInput the trajectory file in .xtc GROMACS format.\n";
	Usage[i++]="\t -pdb <string filename>\n"
			"\t\tInput a PDB file containing the entire system investigated. This input file \n "
			"\t\tis required when computing the SAXS intensity of the solvated protein or of the buffer. From this\n"
			"\t\tfile the topology of the system is extracted.\n";
	Usage[i++]="\t -o <string filename>\n"
			"\t\tOutput a text file containing the spherical averaged SAXS profile from the run.\n";
	Usage[i++]="\t -b <int nstart>\n"
			"\t\tNumber of step from where to begin reading of the system trajectory.This is required when calculating\n"
			"\t\tthe SAXS intensity for a system. \n";
	Usage[i++]="\t -e <int nend>\n"
			"\t\tNumber of step where to end reading the trajectory file. This is required when calculating the\n"
			"\t\tSAXS intensity for a system.\n";
	Usage[i++]="\t -skip <int nskip>\n"
			"\t\tNumber of step  to skip.\n";
	Usage[i++]="\t -solute <string selection>\n"
			"\t\tDefine the residue of the solute to be use in clustering and solute reconstruction. selection is a\n"
			"\t\tstring between double quotes defining the solute residue. If one has a protein as a solute use the\n"
			"\t\tstring ""Protein"" to define as solute all protein residues. This option is required when computing\n"
					"\t\tthe SAXS intensity.\n";
	Usage[i++]="\t -select <string selection>\n"
			"\t\tSelect the system. selection is a doubly quoted string defining the system extent. It must contain\n"
			"\t\tresidues from the solvent and system.\n";

	Usage[i++]="\t -grid <int Nx=128> <int Ny=Nx> <int Nz=Nx>\n"
			"\t\tDefine the direct and reciprocal space grid of the SAXS intensity calculation. It expects at least\n"
			"\t\tNx, or in alternative the three dimensions, Nx, Ny and Nz.\n";
	Usage[i++]="\t -clust <float cutoff [0.0] in Angstroems>\n"
			"\t\tDo clustering by percolation of the solute. This is used to make sure that the solute is at the\n"
			"\t\tcenter of the cell. Doing clustering is expensive, it is always better to have the solute in the\n"
			"\t\tmiddle already at simulation time. If cutoff is zero, the default, the cutoff is chosen according\n"
			"\t\tmiddle to the Lennard-Jones sigma parameter\n";

	Usage[i++]="\t -once // Compute cluster size of detergent heads and ions only once! \n";
	Usage[i++]="\t -dens <string select='R'> <int order=4> <int avg=2>\n"
			"\t\t Compute electron density instead of SAXS. select=[i++] R compute on R-space\n "
			"\t\t Q compute it from Q-space; order is the Lagrangian order; avg is how many bins\n"
			"\t\t are averaged";
	Usage[i++]="\t -saxs <int order=1> <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by linear (order=1) or Langrangian (order >1 ) interpolation. dq is the bin size and qcut is\n"
			"\t\tthe cutoff of the spherical averaged SAXS intensity histogram, printed at the end of each run. dq\n"
			"\t\tand qcut are in Å−1\n";
	Usage[i++]="\t -saxssp <int order=1> <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by cardinal B-spline interpolation. order should be at least 4. Histogram parameters as in -saxs.\n";
	Usage[i++]="\t -obin <string filename>\n"
			"\t\tOutput a binary file containing the 3D SAXS intensity for the system investigated.\n";
	Usage[i++]="\t -in <string filename1> <string filename 2>\n"
			"\t\tTo be used to compute the difference SAXS profile between the SAXS intensity of the solution\n"
			"\t\tfrom file filename1 and that of the buffer filename2\n";
	Usage[i++]="\t -file filename file filename contains commands in free format\n";
	Usage[i++]="\t -padd <string mypadd>\n"
			"\t\tThe trjSaxs default is doing an uniform padding taking the simulation box border density of the\n"
			"\t\tsolvent. This option define alternative padding schemes. The input mypadd can be zero, pbc or a\n"
			"\t\tfilename. In the latter case one has to provide the number density (in Å−3) for each residue of \n"
			"\t\tthe buffer. In the file, place one line for each residue then on each line the name of the residue \n"
			"\t\tfound in the PDB file followed by its number density.\n";
	Usage[i++]="\t -qhist <float dq=0.05> <float qcut=4>\n"
			"\t\tHistogram parameters as in -saxs, but for the difference calculation.\n";

	Usage[i++]="\t -debey <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by Debey formula. Only applicable to finite systems, and very slow. Histogram parameters as in\n"
			"\t\t-saxs.\n";
	Usage[i++]="\t -supcell <float sigma>\n"
			"\t\tDefine the SAXS calculation super-cell. sigma is the multiplicative factor.\n";
	Usage[i++]="\t -direct <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by Debey histogram formula. Only applicable to finite systems, faster than -debey.\n"
			"\t\tHistogram parameters as in -saxs.\n";
	Usage[i++]="\t -nosplout Do not use Hakima spline to print I(q) away from histogram nodes\n";
	Usage[i++]="\t -molw <float molw>\n"
			"\t\tMolecular weight of the solute. molw is the molecular weight in Daltons. Uses the Grishaev et al.\n"
			"\t\tformula to obtain the scaling factor α.\n";

	Usage[i++]="\t -nofluct <int nacc> Run calculation without fluctuations. Averaged on nacc steps\n";
	Usage[i++]="\t -nosolv Run calculation without solvent.\n";
	Usage[i++]="\t -what <sq,[saxs],sans,el> do SQ/SAXS/SANS/ELDENS calculation \n";
	Usage[i++]="\t -lagrange <int order=1> <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by linear (order=1) or Langrangian (order >1 ) interpolation. dq is the bin size and qcut is\n"
			"\t\tthe cutoff of the spherical averaged SAXS intensity histogram, printed at the end of each run. dq\n"
			"\t\tand qcut are in Å−1\n";
	Usage[i++]="\t -spline <int order=1> <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by cardinal B-spline interpolation. order should be at least 4. Histogram parameters as in -saxs.\n";
	Usage[i++]="\t -h // write this help \n";

//	Usage[56]="\t -help // write some on line help \n";
	// Delete inactive commands

	if(ntot <2){
		return ;
	}
	vector<string> vv0;
	string vv;
	if(string(v[1]) =="-file"){

		std::ifstream fin;
		fin.open(v[2],std::ios::in);
		copy(std::istream_iterator<string>(fin),std::istream_iterator<string>(),std::back_inserter(vv0));
	} else {
		for(auto o=1;o<ntot;o++) vv0.push_back(v[o]);
	}
	string key;

	for(auto tmp0: vv0){
		if(tmp0 =="-h" || tmp0 == "--help") {
			for(auto ss: this->getUsage() )
				cout << ss;
			cout << endl;
			exit(1);
		}
		if(tmp0[0] =='-'){
			key.assign(tmp0);
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
	}
	// Eliminate -select with quote problem, which is split by copy

	if(string{v[1]} == "-file"){
		vector<string> sst{"-select","-solute"};

		for(auto myStr: sst){
			if(!inmap[myStr].empty()){
				vector<string> tmp;
				string str;
				for(const auto & op: inmap[myStr]){
					if(op[0]=='-')
						tmp.push_back(op);
					else{
						str+=op+" ";
					}
				}
				inmap[myStr].erase(inmap[myStr].begin(),inmap[myStr].end());
				str.erase(std::remove(str.begin(), str.end(), '"'), str.end());
				tmp.push_back(str);
				inmap[myStr]=tmp;
			}
		}
	}
}

vector<string> trjInput::getUsage(){
	vector<string> use;
	auto it=Usage.begin();
	for(;it!=Usage.end();++it){
		if(!Usage.empty())
			use.push_back(it->second);
	}
	return use;
}
trjInput::~trjInput() {
	// TODO Auto-generated destructor stub
}

} /* namespace trj */
