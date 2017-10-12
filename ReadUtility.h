#ifndef READUTILITY_H
#define READUTILITY_H


#include "GeneralSet.h"
#include"KmerUtility.h"
#include <vector>
#include <string>


using namespace std;


void load_reads(string file,bool rev);

/*
class InfoReadset  {

public:

	InfoReadset() : count(0) {}
	void add(size_t i);
	vector<size_t> get_readset();
	size_t get_count();
	void clear();

private:

	size_t count;
	vector<size_t> readset;

};
*/

#endif
