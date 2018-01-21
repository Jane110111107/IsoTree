// This file was modified from a file named kmerhash.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpacker Software LICENSE.



#include "ReadHash.h"
#include <algorithm>
#include <iostream>
#include <math.h>
#include <time.h>


using namespace std;




vector<int>& ReadHash::operator[](read_int_type read_int) {

	return read_hash[read_int];

}


void ReadHash::get_seeds(vector<int>& seeds) {
	
	read_hash_type_iterator it;
	for (it = read_hash.begin(); it != read_hash.end(); ++it) {
		if (it->second.size() > g_min_seed_coverage && it->second.size() < g_min_seed_coverage+2) {
			seeds.push_back(it->second[0]);
			data_tag[it->second[0]] = -5;
		}
	}	
}

void ReadHash::get_read_hash() {
	
	if (data.empty())
		return;
	time_t start_time = time(NULL);
	size_t data_size = data.size();

	for (int i = 0; i < data_size; ++i) {
		const string& read = data[i];
		read_int_type read_int = get_read_int_type(read);
		read_hash[read_int].push_back(i);

	}
	time_t end_time = time(NULL);
	cout << "Read Hash has been constructed (times: "<< (end_time - start_time) << " s)" << endl;
	read_hash_type_iterator it;
	int total = 0;
	for (it = read_hash.begin(); it != read_hash.end(); ++it)
		++total;
	cout << "There are total of " << total << " read hash" << endl;

	g_min_trans_cov = data_size * g_read_length * 0.05 / total;
	g_min_seed_coverage = data_size / total;

	cout << "seed cov: " << g_min_seed_coverage << " min trans: " << g_min_trans_cov << endl;


}

int ReadHash::unused_sum(read_int_type read_int) {

	int total_unused = 0;
	vector<int>& unused = read_hash[read_int];
	for (int i = 0; i < unused.size(); ++i) {
		if (data_tag[unused[i]] == -1)
			++ total_unused;
	}
	return total_unused;
}


int ReadHash::unused_sum1(read_int_type read_int) {

	int total_unused = 0;
	vector<int>& unused = read_hash[read_int];
	for (int i = 0; i < unused.size(); ++i) {
		if (data_tag[unused[i]] != -2)
			++ total_unused;
	}
	return total_unused;
}



void ReadHash::return_unused(read_int_type read_int, vector<int>& unused_set) {

	unused_set.clear();
	vector<int> read_set = read_hash[read_int];
	for (int i = 0; i < read_set.size(); ++i) {
		if (data_tag[read_set[i]] == -1)
			unused_set.push_back(read_set[i]);
	}
}


