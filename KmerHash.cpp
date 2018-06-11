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

void KmerHash::get_seed(vector<int>& seeds) {

	kmer_hash_type_iterator it;
	vector<kmer_int_type> seed_kmer;

	for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {
		if ( static_cast<int>((it->second).size()) < g_min_seed_coverage || compute_entropy(it->first, kmer_length) < g_min_seed_entropy)
			continue;
		seed_kmer.push_back(it -> first);
	}
	kmer_sorter_by_count_desc_t sorter(*this);
	sort(seed_kmer.begin(), seed_kmer.end(), sorter);
	vector<int> current_seeds;
	for (int i = 0; i < seed_kmer.size(); ++i) {
		current_seeds = kmer_hash[seed_kmer[i]];
		for (int j = 0; j < current_seeds.size(); ++j) {
			if (data_tag[current_seeds[j]] == -4) {
				seeds.push_back(current_seeds[j]);
				data_tag[current_seeds[j]] = -5;
			}
		}
	}
}

void KmerHash::get_left_hash() {


	size_t data_size = data.size();

	cout << "constructing left kmer hash ..." << endl;

	time_t start_time = time(NULL);

	for (size_t i = 0; i < data_size; i++) {

		const string& read = data[i];
		const string& kmer = read.substr(0, kmer_length);
		kmer_int_type kmer_int = kmer_to_int(kmer);
		kmer_hash[kmer_int].push_back(i);
	}

	time_t end_time = time(NULL);
	cout << "Left Kmer Hash has been constructed, total " << kmer_hash.size() << " kmers! (elapsed time: " << (end_time - start_time) << " s)" << endl;

}


void KmerHash::get_right_hash() {


	size_t data_size = data.size();

	cout << "constructing right kmer hash ..." << endl;

	time_t start_time = time(NULL);

	for (size_t i = 0; i < data_size; i++) {

		const string& read = data[i];
		const string& kmer = read.substr(read.length()-kmer_length);
		kmer_int_type kmer_int = kmer_to_int(kmer);
		kmer_hash[kmer_int].push_back(i);
	}

	time_t end_time = time(NULL);
	cout << "Right Kmer Hash has been constructed, total " << kmer_hash.size() << " kmers! (elapsed time: " << (end_time - start_time) << " s)" << endl;

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

	for (it = kmer_hash.begin(); it != kmer_hash.end(); ++it) {

		kmer_int_type kmer = it -> first;
		vector<kmer_pair> candidates; 

		get_forward_candidates_dele(kmer,candidates);

		int max_count = 0;
		int dominant_count = 0;
		for( unsigned int i = 0; i < candidates.size(); ++i) {

			if (candidates[i].second) {
				
				int candidate_count = candidates[i].second;
				if (dominant_count == 0)
					dominant_count = candidate_count;
				else {
					if ( (float) candidate_count/dominant_count < min_ratio_non_error ) {
						kmer_hash_type_iterator candidate = find_kmer(candidates[i].first);
						delete_list.push_back(candidate->first);
						candidate -> second.clear();	
					}
				}

			}

		}
	}
	if (!delete_list.empty()) {

		for (unsigned int i = 0; i < delete_list.size(); ++i) {

			vector<int> dele_reads = kmer_hash[delete_list[i]];
			vector<int>::iterator its;
			for (its = dele_reads.begin(); its != dele_reads.end(); ++its)
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
	vector<int> dele_reads;
	for (pos = kmer_hash.begin(); pos != kmer_hash.end(); ) {

		string kmer = kmer_to_string(pos->first,kmer_length);

		if (static_cast<int>((pos->second).size()) < min_kmer_coverage || compute_entropy(kmer) < min_kmer_entropy) { 
			dele_reads = pos -> second;
			for (int i = 0; i < dele_reads.size(); ++i) {
				data_tag[dele_reads[i]] = -2;
			}
			kmer_hash.erase(pos++);

		} else {
			++pos;
		}
	}
	time_t end_time = time(NULL);
	cout << "Pruning success! (elapsed time: " << (end_time - start_time) << " s)" <<endl;
	cout << " Now, there are total of " << kmer_hash.size() << " kmers." << endl;

	return (true);
}
