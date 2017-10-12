// This file was modified from a file named kmerhash.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpacker Software LICENSE.



#include "KmerHash.h"
#include <algorithm>
#include <iostream>
#include <math.h>
#include <time.h>


using namespace std;


bool cmp_by_value(const pair<int,int>& a, const pair<int, int>& b) {
	return a.second > b.second;
}

//functions in class KmerHash
KmerHash::KmerHash (int kmer_length) {

	this -> kmer_length = kmer_length;
	

}

//template<typename T>
vector<int>& KmerHash::operator[](kmer_int_type kmer) {

	return kmer_hash[kmer];

}


//template<typename T>
size_t KmerHash::get_size() {

	return kmer_hash.size();

}


//template<typename T>
bool KmerHash::empty() {

	return kmer_hash.empty();

}


//template<typename T>
bool KmerHash::reuse(const kmer_int_type kmer) {

	kmer_hash_type_iterator it = find_kmer(kmer);
	return (it != kmer_hash.end());
	
}


//template<typename T>
bool KmerHash::exists(kmer_int_type kmer) {

	return (kmer_abundance(kmer) > 0);
	
}


//template<typename T>
bool KmerHash::exists(const string& kmer) {

	kmer_int_type it = kmer_to_int(kmer);
	return (exists(it));
	
}


//template<typename T>
size_t KmerHash::kmer_abundance(kmer_int_type kmer) {

	kmer_hash_type_iterator it = find_kmer(kmer);

	if (it != kmer_hash.end())
		return ((it -> second).size());
	else
		return 0;

}


//template<typename T>
size_t KmerHash::kmer_abundance(const string& kmer) {

	kmer_int_type it = kmer_to_int(kmer);
	return (kmer_abundance(it));

}


void KmerHash::get_hash(vector<int>& seeds) {

	size_t data_size = data.size();

	cout << "constructing kmer hash ..." << endl;

	time_t start_time = time(NULL);
	
	if (data.empty())
		return;

	for (size_t i = 0; i < data_size; i++) {

		const string& read = data[i];

		if ((int)read.length() < kmer_length) 
			continue;
		if (contains_non_gatc(read)) {
			data_tag[i] = -2;
			int mate_id;
			if (i >= max_read_id)
				mate_id = i - max_read_id;
			else
				mate_id = i + max_read_id;
			continue;
		}
		if (data_tag[i] == -2)
			continue;
		bool is_seed = true;
		for (size_t j = 0; j <= read.length()-kmer_length; j++) {

			const string& kmer = read.substr(j, kmer_length);			
			kmer_int_type kmer_int = kmer_to_int(kmer);

			if (g_is_paired_end && g_double_stranded_mode) {
				kmer_int_type kmer_rev = revcomp_val(kmer_int,kmer_length);
				if (exists(kmer_rev)) {
					kmer_hash[kmer_rev].push_back(i);
					is_seed = false;
				}
				else {
					vector<int> info;
                    			info.push_back(i);
			        	kmer_hash[kmer_rev] = info;
					if (kmer_hash.size() % 1000000 == 0)
						cout << "constructing " << kmer_hash.size() << " kmers..." << endl;

			        //if (j == 0 && tag_first && i < max_read_id) {				    
						//seeds.push_back(i);				   
						//tag_first = false;			     
					//} 
				}

			}

			if (exists(kmer_int)) {

				kmer_hash[kmer_int].push_back(i);
				is_seed = false;
				continue;

			}

			vector<int> info;
            info.push_back(i);
			kmer_hash[kmer_int] = info;
			if (kmer_hash.size() % 1000000 == 0)
				cout << "constructing " << kmer_hash.size() << " kmers..." << endl;

			/*if (j == 0 && tag_first) {
				if (g_is_paired_end && i < max_read_id) {
					seeds.push_back(i);
					tag_first = false;
				}
				if (!g_is_paired_end) {
					seeds.push_back(i);
					tag_first = false;
				}
			} 
			*/
		}

		if(is_seed) {
			//if (g_is_paired_end && i < max_read_id)
				//seeds.push_back(i);
			//if (!g_is_paired_end)
				seeds.push_back(i);
		}
	}
/*
	cout << "get sorted seeds..." << endl;
	vector<pair<int, int> > seed_sum;
	for (int i = 0; i < seeds.size(); i++) {
		int x = seeds[i];
		int ss = 0;
		for (int j = 0; j < g_read_length-g_kmer_length; j++) {
			string k_str = data[x].substr(j, g_kmer_length);
			kmer_int_type kmer_int = kmer_to_int(k_str);
			ss = ss + kmer_hash[kmer_int].size();
		}
		pair<int, int> s_t;
		s_t.first = x;
		s_t.second = ss;
		seed_sum.push_back(s_t);
	}
cout << seeds.size() << endl;
cout << seed_sum.size() << endl;
cout << "seed sort" << endl;

sort(seed_sum.begin(), seed_sum.end(), cmp_by_value);
for (int i = 0; i < seed_sum.size(); i++)
	seeds[i] = seed_sum[i].first;
*/
	time_t end_time = time(NULL);
	cout << "Kmer Hash has been constructed, total " << kmer_hash.size() << " kmers! (elapsed time: " << (end_time - start_time) << " s)" << endl;

}


void KmerHash::get_hash_from_graphdata(vector<string>& mdata) {

	size_t mdata_size = mdata.size();
	
	if (mdata.empty())
		return;

	for (size_t i = 0; i < mdata_size; i++) {

		const string& read = mdata[i];

		if ((int)read.length() < kmer_length) 
			continue;

		for (size_t j = 0; j <= read.length()-kmer_length; j++) {

			const string& kmer = read.substr(j, kmer_length);

			if (contains_non_gatc(kmer))

				continue;

			kmer_int_type kmer_int = kmer_to_int(kmer);

			if (exists(kmer_int)) {

				kmer_hash[kmer_int].push_back(i);
				continue;

			}

			vector<int> info;
            info.push_back(i);
			kmer_hash[kmer_int] = info;

		}

	}

}


//template<typename T>
bool KmerHash::remove(kmer_int_type kmer) {

	kmer_hash_type_iterator it = find_kmer(kmer);

	if (it != kmer_hash.end()) {

		kmer_hash.erase(it);
		return (true);

	}

	return (false);

}



void KmerHash::get_forward_candidates_dele(kmer_int_type seed_kmer, vector<kmer_pair>& candidates) {

	candidates.clear();
	kmer_int_type forward_prefix = (seed_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;
	for (kmer_int_type i = 0; i < 4; i++) {
		kmer_pair candidate;
		candidate.first = forward_prefix | i;
		candidate.second = kmer_abundance(candidate.first);
		if (candidate.second)
			candidates.push_back(candidate);
	}

}


bool KmerHash::delete_errous_kmer(float min_ratio_non_error) {

	time_t start_time = time(NULL);
	cout << " Remove erroneous kmers..." << endl;
	kmer_hash_type_iterator it;
	vector<kmer_int_type> delete_list;

	for (it = kmer_hash.begin(); it != kmer_hash.end(); it++) {

		kmer_int_type kmer = it -> first;
		vector<kmer_pair> candidates; 

		get_forward_candidates_dele(kmer,candidates);

		int max_count = 0;

		for( unsigned int i = 0; i < candidates.size(); i++) {

			if (candidates[i].second) {

				int candidate_count = candidates[i].second;

				if (max_count < candidate_count)

					max_count = candidate_count;

			}

		}

		for( unsigned int i = 0; i < candidates.size(); i++) {

			if (candidates[i].second) {

				int candidate_count = candidates[i].second;

				if ( (float) candidate_count/max_count < min_ratio_non_error) {

					kmer_hash_type_iterator candidate = find_kmer(candidates[i].first);
					delete_list.push_back(candidate->first);
					candidate -> second.clear();
				}

			}

		}

	}


	if (!delete_list.empty()) {

		for (unsigned int i = 0; i < delete_list.size(); i++) {

			vector<int> dele_reads = kmer_hash[delete_list[i]];
			vector<int>::iterator its;
			for (its = dele_reads.begin(); its != dele_reads.end(); its++)
				data_tag[*its] = -2;

			remove(delete_list[i]);
		}

		time_t end_time = time(NULL);
		cout << "Remove success! (elapsed time: " << (end_time - start_time) << " s)" <<endl;
		cout << " Now, there are total of " << kmer_hash.size() << " kmers." << endl;
		return (true);

	}

	return (false);

}


bool KmerHash::prune_hash(int min_kmer_coverage, float min_kmer_entropy) {

	if (min_kmer_coverage == 1 && min_kmer_entropy < 1e-6)
		return (true);

	time_t start_time = time(NULL);

	cout << "Pruning kmer hash ..." << endl;

	kmer_hash_type_iterator pos;
	vector<kmer_int_type> delete_list;

	for (pos = kmer_hash.begin(); pos != kmer_hash.end(); pos++) {

		string kmer = kmer_to_string(pos->first,kmer_length);

		if (static_cast<int>((pos->second).size()) < min_kmer_coverage || compute_entropy(kmer) < min_kmer_entropy) 
			delete_list.push_back(pos->first);

	}

	time_t end_time = time(NULL);

	if (!delete_list.empty()) {

		for (unsigned int i = 0; i < delete_list.size(); i++) {

			vector<int> dele_reads = kmer_hash[delete_list[i]];
			vector<int>::iterator its;
			for (its = dele_reads.begin(); its != dele_reads.end(); its++)
				data_tag[*its] = -2;

			remove(delete_list[i]);

		}
		time_t end_time = time(NULL);
		cout << "Pruning success! (elapsed time: " << (end_time - start_time) << " s)" <<endl;
		cout << " Now, there are total of " << kmer_hash.size() << " kmers." << endl;
		return (true);

	}

	return (false);

}
