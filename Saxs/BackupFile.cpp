/*
 * BackupFile.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: marchi
 */

#include "BackupFile.h"
auto checkFile=[](string s)-> bool {
	ifstream f(s);
	if(f.good()) {
		f.close();return true;
	} else {
		f.close();return false;}
};
BackupFile::BackupFile(string ss) {
	ifstream f(ss.c_str());
	if(f.good()){
		// find if backup exists
		int n{0};
		FilenameBackup s0(ss);
		string s1{s0(n)};
		while(checkFile(s1)){
			n++;
			s1=s0(n);
		}
		Name=s1;
		ifstream source(ss.c_str(), ios::binary);
	    ofstream dest(Name.c_str(), ios::binary);
	    cout << "\n<---- Backing up binary output file to " + Name + " ---->"<<endl;
	    dest << source.rdbuf();
	    source.close();
	    dest.close();

	}
	f.close();
}

BackupFile::~BackupFile() {
	// TODO Auto-generated destructor stubw
}

