//
// win2unix.cpp
//
// Author: Jennifer Liddle (js10)
//
// $Id: win2unix.cpp 1241 2010-07-27 16:14:37Z scp $
//
//
#include "win2unix.h"
#include <iostream>

using namespace std;

string filepath;

string findAndReplace(string target, string find, string replace)
{
    size_t pos = target.find(find);
    while (pos != string::npos) {
        target.replace(pos,find.length(),replace);
        // unlikely we'll exchange strings which are substrings of one
        // another. just in case, start next search at the next
        // char (pos + 1) + offset (difference in string lengths)
        pos = target.find(find, pos + 1 + replace.length() - find.length());
    }
    return target;
}


const char * win2unix(const char *f)
{
	filepath = f;
	filepath = findAndReplace(filepath,"\\","/");
	filepath = findAndReplace(filepath,"//","/");
	filepath = findAndReplace(filepath,"evs-users3","evs-illumina");
	filepath = findAndReplace(filepath,"evs-illumina","nfs");
	filepath = findAndReplace(filepath,"fastnfs/illumina","nfs/new_illumina");
	filepath = findAndReplace(filepath,"geno1","geno01");
	filepath = findAndReplace(filepath,"geno2","geno02");
	filepath = findAndReplace(filepath,"geno3","geno03");
	filepath = findAndReplace(filepath,"geno4","geno04");
	filepath = findAndReplace(filepath,"geno5","geno05");
	filepath = findAndReplace(filepath,"geno6","geno06");
	filepath = findAndReplace(filepath,"geno7","geno07");
	filepath = findAndReplace(filepath,"geno8","geno08");
	filepath = findAndReplace(filepath,"geno9","geno09");
	return filepath.c_str();
}

string win2unix(string f)
{
	filepath = win2unix(f.c_str());
	return filepath;
}


