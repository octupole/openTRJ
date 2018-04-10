/*
 * SaxsPadding.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: marchi
 */

#include "SaxsPadding.h"

void SaxsPadding::operator()(vector<string> & vstr,vector<double> vint){
	List=vstr;
	iList=vint;
}
void SaxsPadding::CheckIt(vector<string> & Ref, vector<string> & Solu){
	if(!List.size()) {Have_CheckedIt=true;return;}
	sort(Ref.begin(),Ref.end());
	sort(Solu.begin(),Solu.end());
	vector<string> my{List},out,out0;
	sort(my.begin(),my.end());
	try{
		set_intersection(my.begin(),my.end(),Ref.begin(),Ref.end(),back_inserter(out));
//		if(List.size() != out.size()) throw string("\nAll padding species must be found in the SAXS selection.\n");
		out.clear();
		set_intersection(my.begin(),my.end(),Solu.begin(),Solu.end(),back_inserter(out));
		if(out.size()) throw string("\nSome padding species are in the solute definition -solute. This cannot be.\n");

		// Add the additional selection types with iList o
		out.clear();
		out0.clear();
		set_difference(Ref.begin(),Ref.end(),Solu.begin(),Solu.end(),back_inserter(out0));
		set_difference(out0.begin(),out0.end(),my.begin(),my.end(),back_inserter(out));
		for(auto op: out) {
			List.push_back(op);
			iList.push_back(0.0);
		}

	}catch(const string & s){cout << s <<endl;Finale::Finalize::Final();}
	Have_CheckedIt=true;
}
void SaxsPadding::setMapResidue(const map<string,vector<string> > & mapp){
	try{
		if(!Have_CheckedIt) throw string("Must Check the padding data before usage!!");
	}catch(const string & s){cout << s <<endl;Finale::Finalize::Final();}
	for(auto o=0;o<List.size();o++){
		auto y=mapp.find(List[o]);
		map<string,double> tmp;
		for(auto opy: y->second){
			tmp[opy]+=iList[o];
		}
		for(auto it=tmp.begin();it != tmp.end();++it){
			MapResidue[it->first]+=it->second;
		}
	}
 	cout << std::fixed << std::setw(6)<< "\n"+string(80,'-')+" \n\n            Padding of the repeated cells"
 			<< "\n\n";
 	cout << std::showpoint << std::fixed << std::left;
 	cout << std::setw(20) << "       Atom Name"<< std::setw(20) << "Cell Number Density [1/nm^{3}]\n" << endl;

 	size_t AtomTotal=0;
	for(auto it=MapResidue.begin();it != MapResidue.end();++it){
 		cout << std::left << std::fixed<< "         " <<std::setw(20) <<
 				it->first << "      " << std::setw(24) << it->second
 				<< endl;
	}
	cout << std::fixed << std::setw(6)<< "\n"+string(80,'-') +"\n";
}
SaxsPadding::~SaxsPadding() {
	// TODO Auto-generated destructor stub
}

