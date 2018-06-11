// This file was modified from a file named sequenceUtil.hpp in Inchworm modules of Trintiy.
// The original copyright info is Copyright (c) 2010, The Broad Institute, Inc.
// Distributed under the  Distributed under Trinity Software LICENSE.



#ifndef KMERUTILTIY_H
#define KMERUTILITY_H



#include <iostream>
#include <cstring>
#include <cctype>
#include <map>
#include <string>
#include <vector>


using namespace std;


typedef unsigned long long kmer_int_type;


//check a kmer wether contains characters that not G/A/T/C
bool contains_non_gatc(const string& kmer);


char int_to_base(int baseval); // 0 1 2 3 => G A T C
int base_to_int(char uncleotide); //(GATC) = {0 1 2 3}, others = -1


//convert a kmer in string format into a 64 interger
kmer_int_type kmer_to_int(const string& kmer);//must be less than 32 bases


//covert a kmer of a 64 interger into a string
string kmer_to_string(kmer_int_type kmer, unsigned int kmer_length);


//compute entropy
float compute_entropy(const string& kmer);
float compute_entropy(kmer_int_type kmer, unsigned int kmer_length);


//return the complementary of a kmer in reverse order
string revcomp(const string& kmer);
kmer_int_type revcomp_val(kmer_int_type kmer, unsigned int kmer_length);


//store a kmer as itself or its reverse complement
kmer_int_type get_DS_kmer_val(kmer_int_type kmer, unsigned int kmer_length);


#endif
