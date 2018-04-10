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
	inmap["-prof"]=in;
	inmap["-pVct"]=in;
	inmap["-pdx1"]=in;
	inmap["-pdx2"]=in;
	inmap["-cut"]=in;
	inmap["-rcut"]=in;
	inmap["-cutX"]=in;
	inmap["-grid"]=in;
	inmap["-chg"]=in;
	inmap["-rho"]=in;
//	inmap["-help"]=in;
	inmap["-irl"]=in;
	inmap["-irdb"]=in;
	inmap["-rd"]=in;
	inmap["-ic"]=in;
	inmap["-clust"]=in;
	inmap["-diffk"]=in;
	inmap["-diffq"]=in;
	inmap["-diffm"]=in;
	inmap["-Voro"]=in;
	inmap["-rin"]=in;
	inmap["-filter"]=in;
	inmap["-bsp"]=in;
	inmap["-in"]=in;
	inmap["-in2"]=in;
	inmap["-e0"]=in;
	inmap["-file"]=in;
	inmap["-neq"]=in;
	inmap["-dip"]=in;
	inmap["-lagrange"]=in;
	inmap["-spline"]=in;
	inmap["-saxs"]=in;
	inmap["-saxssp"]=in;
	inmap["-what"]=in;
	inmap["-lowq"]=in;
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

	map<string,vector<string> >::iterator it=inmap.begin();
	for(int n=0;it!=inmap.end();++it,n++){
		Usage[n]=" ";
	}
	Usage[0]="\t -dcd <string filename>\n"
			"\t\tInput the trajectory file in NAMD format. This is experimental and \n"
			"\t\thas not yet been thoroughly tested.\n";
	Usage[1]="\t -xtc <string filename> \n"
			"\t\tInput the trajectory file in .xtc GROMACS format.\n";
	Usage[2]="\t -pdb <string filename>\n"
			"\t\tInput a PDB file containing the entire system investigated. This input file \n "
			"\t\tis required when computing the SAXS intensity of the solvated protein or of the buffer. From this\n"
			"\t\tfile the topology of the system is extracted.\n";
	Usage[3]="\t -o <string filename>\n"
			"\t\tOutput a text file containing the spherical averaged SAXS profile from the run.\n";
	Usage[4]="\t -b <int nstart>\n"
			"\t\tNumber of step from where to begin reading of the system trajectory.This is required when calculating\n"
			"\t\tthe SAXS intensity for a system. \n";
	Usage[5]="\t -e <int nend>\n"
			"\t\tNumber of step where to end reading the trajectory file. This is required when calculating the\n"
			"\t\tSAXS intensity for a system.\n";
	Usage[6]="\t -skip <int nskip>\n"
			"\t\tNumber of step  to skip.\n";

	Usage[7]="\t -nohyd // Do not include hydrogens in Micelle calculation [default] \n";
	Usage[8]="\t -hyd // Include hydrogens in the Micelle calculation \n";

	Usage[9]="\t -nodel // Do not delete temporary files in MPI runs \n";
	Usage[10]="\t -del // Delete temporary files in MPI runs [default]\n";
	Usage[11]="\t -solute <string selection>\n"
			"\t\tDefine the residue of the solute to be use in clustering and solute reconstruction. selection is a\n"
			"\t\tstring between double quotes defining the solute residue. If one has a protein as a solute use the\n"
			"\t\tstring ""Protein"" to define as solute all protein residues. This option is required when computing\n"
					"\t\tthe SAXS intensity.\n";
	Usage[12]="\t -select <string selection>\n"
			"\t\tSelect the system. selection is a doubly quoted string defining the system extent. It must contain\n"
			"\t\tresidues from the solvent and system.\n";
	Usage[13]="\t -prof // Compute the radial and 2D density profiles.\n";
	Usage[14]="\t -cut cutoff // General cutoff. For the radial profile cutoff is the spherical cuoff.\n"
			"\t\tFor the 2D profile cutoff is the on plane cutoff. Default value is 40 A\n";
	Usage[15]="\t -cutX cutoff // Out of plane cutoff for 2D profiling. Default value is 5 A.\n";
	Usage[16]="\t -pVct (string) Axes// Axes can take values as \"X\", \"Y\" or \"Z\". \n  "
			"\t\tThis is the axis perpendicular to the plane to be profiled to output. Default is \"X\".\n";
	Usage[17]="\t -pdx1 dx1 // Bin size of the radial histograms.\n";
	Usage[18]="\t -pdx2 dx2 // Bin size of the 2D histograms.\n";
	Usage[19]="\t -chg filename // Input charges for the selected residues.\n";
	Usage[20]="\t -rho[order,dq,qcut] do Electrostatics calculation with Lagrangian interpolation. If order = 0 use a tridiagonal approximation\n";
	Usage[21]="\t -grid <int Nx=128> <int Ny=Nx> <int Nz=Nx>\n"
			"\t\tDefine the direct and reciprocal space grid of the SAXS intensity calculation. It expects at least\n"
			"\t\tNx, or in alternative the three dimensions, Nx, Ny and Nz.\n";
	Usage[22]="\t -irl   // Input radius list to be used with the -w option. Default values:\n"
			"\t\t Tsai et al., J Mol Biol 290:253 (1999) \n";
	Usage[23]="\t -idb   // Input residue database for radius to be used with the -w option \n";
	Usage[24]="\t -rd // Use atomic radii in the Micelle calculation \n";
	Usage[25]="\t -clust <float cutoff [0.0] in Angstroems>\n"
			"\t\tDo clustering by percolation of the solute. This is used to make sure that the solute is at the\n"
			"\t\tcenter of the cell. Doing clustering is expensive, it is always better to have the solute in the\n"
			"\t\tmiddle already at simulation time. If cutoff is zero, the default, the cutoff is chosen according\n"
			"\t\tmiddle to the Lennard-Jones sigma parameter\n";

	Usage[26]="\t -rcut cutoff // Cutoff for Linked Lists. Default value is 35 A.\n";
	Usage[27]="\t -once // Compute cluster size of detergent heads and ions only once! \n";
	Usage[28]="\t -rcm file.rcm\n";
	Usage[29]="\t -Voro Use Voronoi to compute neighbors. More expensive, but sure method\n";
	Usage[30]="\t -rin Inner cutoff in A for neighbor list. By default 30, 15 A in case -Voro \n";
	Usage[31]="\t -diffk  do atomic diffusion \n";
	Usage[32]="\t -diffq do rotational diffusion for detergent molecules. The three Wigner arguments in a string format are in input\n";
	Usage[33]="\t -diffm do rotational and translational diffusion for the micelles. Same arguments as for -diffq\n";
	Usage[34]="\t -saxs <int order=1> <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by linear (order=1) or Langrangian (order >1 ) interpolation. dq is the bin size and qcut is\n"
			"\t\tthe cutoff of the spherical averaged SAXS intensity histogram, printed at the end of each run. dq\n"
			"\t\tand qcut are in Å−1\n";
	Usage[35]="\t -saxssp <int order=1> <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by cardinal B-spline interpolation. order should be at least 4. Histogram parameters as in -saxs.\n";
	Usage[36]="\t -ic fileContrast \n";
	Usage[37]="\t -lowq nR [3] Do 3D padding, by adding nR cells on each direction of the saxs cell\n";
	Usage[38]="\t -filter Sigma // use filter with a given gaussian width \n";
	Usage[39]="\t -bsp alpha // use B spline alpha is the Ewald sum alpha parameter in A^{-1}\n";
	Usage[40]="\t -xtc2 file2.xtc\n";
	Usage[41]="\t -obin <string filename>\n"
			"\t\tOutput a binary file containing the 3D SAXS intensity for the system investigated.\n";
	Usage[42]="\t -in <string filename1> <string filename 2>\n"
			"\t\tTo be used to compute the difference SAXS profile between the SAXS intensity of the solution\n"
			"\t\tfrom file filename1 and that of the buffer filename2\n";
	Usage[43]="\t -in2 Read a binary input of I(q) for the reference system \n";
	Usage[44]="\t -e0 Applied Electric field in units of V/nm a real vector\n";
	Usage[45]="\t -file filename file filename contains commands in free format\n";
	Usage[46]="\t -neq //Does non equilibrium calculations\n";
	Usage[47]="\t -padd <string mypadd>\n"
			"\t\tThe trjSaxs default is doing an uniform padding taking the simulation box border density of the\n"
			"\t\tsolvent. This option define alternative padding schemes. The input mypadd can be zero, pbc or a\n"
			"\t\tfilename. In the latter case one has to provide the number density (in Å−3) for each residue of \n"
			"\t\tthe buffer. In the file, place one line for each residue then on each line the name of the residue \n"
			"\t\tfound in the PDB file followed by its number density.\n";
	Usage[48]="\t -qhist <float dq=0.05> <float qcut=4>\n"
			"\t\tHistogram parameters as in -saxs, but for the difference calculation.\n";

	Usage[49]="\t -dip Dielectric constant calculation \n";
	Usage[50]="\t -debey <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by Debey formula. Only applicable to finite systems, and very slow. Histogram parameters as in\n"
			"\t\t-saxs.\n";
	Usage[51]="\t -supcell <float sigma>\n"
			"\t\tDefine the SAXS calculation super-cell. sigma is the multiplicative factor.\n";
	Usage[52]="\t -direct <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by Debey histogram formula. Only applicable to finite systems, faster than -debey.\n"
			"\t\tHistogram parameters as in -saxs.\n";
	Usage[53]="\t -nosplout Do not use Hakima spline to print I(q) away from histogram nodes\n";
	Usage[54]="\t -molw <float molw>\n"
			"\t\tMolecular weight of the solute. molw is the molecular weight in Daltons. Uses the Grishaev et al.\n"
			"\t\tformula to obtain the scaling factor α.\n";

	Usage[55]="\t -nofluct <int nacc> Run calculation without fluctuations. Averaged on nacc steps\n";
	Usage[56]="\t -nosolv Run calculation without solvent.\n";
	Usage[57]="\t -what <sq,[saxs],sans> do SQ/SAXS/SANS calculation \n";
	Usage[58]="\t -lagrange <int order=1> <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by linear (order=1) or Langrangian (order >1 ) interpolation. dq is the bin size and qcut is\n"
			"\t\tthe cutoff of the spherical averaged SAXS intensity histogram, printed at the end of each run. dq\n"
			"\t\tand qcut are in Å−1\n";
	Usage[59]="\t -spline <int order=1> <float dq=0.05> <float qcut=4>\n"
			"\t\tSAXS by cardinal B-spline interpolation. order should be at least 4. Histogram parameters as in -saxs.\n";

//	Usage[56]="\t -help // write some on line help \n";
	// Delete inactive commands
	ClrU({7,8,9,10,13,14,15,16,17,18,19,20,22,23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,43,44,46,49});
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
