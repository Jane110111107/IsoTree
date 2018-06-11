// This file was modified from a file named kmerhash.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpacker Software LICENSE.



#include "ReadHash.h"
#include <algorithm>
#include <iostream>
#include <math.h>
#include <time.h>


using namespace std;




int& ReadHash::operator[](read_int_type read_int) {

	return read_hash[read_int];

}


void ReadHash::get_read_hash() {
	
	if (data.empty())
		return;
	time_t start_time = time(NULL);
	size_t data_size = data.size();

	const int vec_len = g_read_length / 32;
	read_int_type read_int(vec_len+1);
	
	int max_read_id = data_size / 2;
	int total = 0;	
	int hash_id;
	vector<string> identified;
	vector<string>(data).swap(identified);
	identified.clear();
	vector<vector<int> > data_index;

	for (int i = 0; i < data_size; ++i) {
		const string& read = data[i];
		get_read_int_type(read,read_int, vec_len);
		hash_id = find_read(read_int);
		if (hash_id > 0){
			data_index[hash_id].push_back(i);
		} else {
			read_hash[read_int] = total;
			vector<int> xxx(1, i);
			data_index.push_back(xxx);
			hash_id = total;
			identified.push_back(read);
			total = total + 1;
		}
		data_tag[i] = hash_id;
		if (i == max_read_id)
			g_divide_right_pos = total;
	}

	const string& read1 = data[max_read_id];
	get_read_int_type(read1,read_int, vec_len);
	g_divide_left_pos = find_read(read_int);
	
	vector<vector<int> >(data_index).swap(data_pair);
	data_pair.clear();
	int mate_id = 0;
	vector<int> pair_vec;
	for (int i = 0; i < data_index.size(); ++i) {
		pair_vec.clear();
		for (int j = 0; j < data_index[i].size(); ++j) {
			if (data_index[i][j] >= max_read_id) {
				mate_id = data_index[i][j] - max_read_id;
			} else {
				mate_id = data_index[i][j] + max_read_id;
			}
			pair_vec.push_back(data_tag[mate_id]);
		}
		data_pair.push_back(pair_vec);
	}

	vector<string>(identified).swap(data);
	vector<string>().swap(identified);
	vector<vector<int> >().swap(data_index);
	vector<vector<int> >(data_pair).swap(data_pair);
	const int data_size1 = data.size();
	vector<int> updata_tag(data_size1, -5);
	vector<int>(updata_tag).swap(data_tag);
	vector<int>().swap(updata_tag);
	vector<int>(data_tag).swap(data_cov);
	for (int i = 0; i < data_pair.size(); ++i) 
		data_cov[i] = data_pair[i].size();

	time_t end_time = time(NULL);
	cout << "There all total of " << total << " unique read sequences (elapsed time: "<< (end_time - start_time) << " s)" << endl;

	g_min_trans_cov = data_size * g_read_length * 0.05 / total;
	//cout << "seed cov: " << g_min_seed_coverage << " min trans: " << g_min_trans_cov << endl;
}


