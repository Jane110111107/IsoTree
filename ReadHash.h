// This file was modified from a file named kmerhash.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpacker Software LICENSE.


#ifndef READHASH_H
#define READHASH_H


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
typedef  boost::unordered_map<read_int_type, vector<int> > read_hash_type;
typedef  boost::unordered_map<read_int_type, vector<int> >::iterator read_hash_type_iterator;
typedef  boost::unordered_map<read_int_type, vector<int> >::const_iterator read_hash_type_const_iterator;



class ReadHash {
/*
private:
	class ReadHashSort {
	public:
		ReadHashSort(ReadHash& x) : read_hash_sort(x) {};
		bool operator() (const read_int_type& i, const read_int_type& j) {
			return ( read_hash_sort[i].size() > read_hash_sort[j].size() );
		}
	private:
		ReadHash& read_hash_sort;
	};
*/
private:
	read_hash_type read_hash;
public:

	ReadHash () { }
	vector<int> & operator[](read_int_type read_int);
	void get_read_hash();
	void get_seeds(vector<int>& seeds);
	int unused_sum(read_int_type read_int);
	int unused_sum1(read_int_type read_int);
	void return_unused(read_int_type read_int, vector<int>& unused_set);

};



#endif
