// This file was modified from a file named kmerhash.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpacker Software LICENSE.


#ifndef KMERHASH_H
#define KMERHASH_H


#include "KmerUtility.h"
#include "ReadUtility.h"
#include "GeneralSet.h"
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <boost/unordered_map.hpp>


using namespace std;


typedef pair<kmer_int_type, size_t> kmer_pair;
typedef pair<int,int> pair_t;


bool cmp_by_value(const pair<int,int>& a, const pair<int, int>& b);


class KmerHash {

private:

	typedef  boost::unordered_map<kmer_int_type, vector<int> > kmer_hash_type;
	typedef  boost::unordered_map<kmer_int_type, vector<int> >::iterator kmer_hash_type_iterator;
	typedef  boost::unordered_map<kmer_int_type, vector<int> >::const_iterator kmer_hash_type_const_iterator;

	class kmer_sorter_by_count_desc_t {
	public:
		kmer_sorter_by_count_desc_t(KmerHash& k) : kmer_hash_(k) {};
		bool operator() (const kmer_int_type& i, const kmer_int_type& j) {
			return ( (kmer_hash_[i].size() > kmer_hash_[j].size()) || (kmer_hash_[i].size() == kmer_hash_[j].size() && i > j) );
		}
	private:
		KmerHash& kmer_hash_;
	};

	kmer_hash_type kmer_hash;
	int kmer_length;

public:

	KmerHash () { }
	KmerHash (int kmer_length);
	vector<int> & operator[](kmer_int_type kmer);
	size_t get_size();
	bool empty(); 

	kmer_hash_type_iterator find_kmer(kmer_int_type kmer){
		//if (g_double_stranded_mode)
			//kmer = get_DS_kmer_val(kmer, kmer_length);
		return kmer_hash.find(kmer);
	}

	bool reuse(const kmer_int_type kmer);
	bool exists(const kmer_int_type kmer); //weather find the kmer in hash table
	bool exists(const string& kmer); 
	size_t kmer_abundance(kmer_int_type kmer);
	size_t kmer_abundance(const string& kmer);
	void get_seed(vector<int>& seeds);
	void get_hash(vector<int>& seeds);
	void get_hash();
	void get_hash_from_graphdata(vector<string>& mdata);
	bool remove (kmer_int_type kmer);
	bool delete_errous_kmer(float min_ratio_non_error);
	bool prune_hash(int min_kmer_coverage, float min_kmer_entropy);
	void get_forward_candidates_dele(kmer_int_type seed_kmer, vector<kmer_pair>& candidates);
};



#endif
