/*
 * BackupFile.h
 *
 *  Created on: Apr 8, 2016
 *      Author: marchi
 */

#ifndef SRC_BACKUPFILE_H_
#define SRC_BACKUPFILE_H_
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
/** \brief A functor constructing a filename of a GROMACS backup file
 *
 */

struct FilenameBackup{
	string myName;
	FilenameBackup(string s): myName(s){};
	string operator()(int n){
		stringstream ss;
		ss<<n;
		return "#"+myName+"."+ss.str()+"#";
	}
};
/*! \brief Make backups of files following the Gromacs scheme
 *
 *
 */
class BackupFile {
	string Name; //!< Original name of the file
public:
	BackupFile()=delete;
/**
 *
 * @param y is the name of the file to be backed up if it exists already
 */
	BackupFile(string y);
	virtual ~BackupFile();
};

#endif /* SRC_BACKUPFILE_H_ */
