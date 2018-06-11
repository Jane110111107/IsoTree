#ifndef READUTILITY_H
#define READUTILITY_H


#include "GeneralSet.h"
#include"KmerUtility.h"
#include <vector>
#include <string>


using namespace std;


typedef vector<unsigned long long> read_int_type;


void load_reads(string file, vector<string>& input_data, bool rev);
void delete_error_reads(vector<string>& input_data);
int count_reads_file(string file);
read_int_type get_read_int_type(const string& read);
void get_read_int_type(const string& read, read_int_type& read_int, const int vec_len);

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
