// This file was modified from a file named sequenceUtil.cpp in Inchworm modules of Trintiy.
// The original copyright info is Copyright (c) 2010, The Broad Institute, Inc.
// Distributed under the  Distributed under Trinity Software LICENSE.



#include"KmerUtility.h"
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<sstream>

using namespace std;

char _int_to_base [4] = {'G', 'A', 'T', 'C'};
unsigned char _base_to_int [256] = {
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //   0
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  20
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  40
		255, 255, 255, 255, 255,   1, 255,   3, 255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, //  60
		255, 255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   1, 255,   3, //  80
		255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   2, 255, 255, 255, // 100
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 120
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 140
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 180
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 200
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 220
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255                      // 240
};


bool contains_non_gatc (const string& kmer) {

	for(unsigned int i =0; i < kmer.size(); i++) {

		unsigned char c = kmer[i];

		if(_base_to_int[c] > 3)
			return(true);

	}

	return(false);
}


char int_to_base (int baseval) {

	if (baseval < 0 || baseval > 3) {
		cout << "Error, baseval out of range 0-3" << endl;
		exit(1);
	}

	return (_int_to_base[baseval]);
}


int base_to_int (char uncleotide) {

	switch (uncleotide) {

	case 'G':
	case 'g':
		return (0);

	case 'A':
	case 'a':
		return (1);

	case 'T':
	case 't':
		return (2);

	case 'C':
	case 'c':
		return (3);

	default:
		return (-1);

	}

}



kmer_int_type kmer_to_int (const string& kmer) {

	if (kmer.length() > 32) {
	cout <<  "Error, kmer length exceeds 32" << endl;
	exit(1);
  }

	kmer_int_type kmer_val = 0;

	for (unsigned int i =0; i < kmer.length(); i++) {

		unsigned char c = kmer[i];
		int val = _base_to_int[c];
		kmer_val = kmer_val << 2;
		kmer_val |= val;

	}

	return (kmer_val);

}


string kmer_to_string (kmer_int_type kmer, unsigned int kmer_length) {

	string kmerstring(kmer_length, ' ');

	for (unsigned int i =1; i <= kmer_length; i++){

		int base_num = kmer & 3; 
		kmerstring[kmer_length-i] = _int_to_base [base_num];
		kmer = kmer >> 2;
	
	}

	return (kmerstring);

}


float compute_entropy (const string& kmer) {

	map<char,int> char_map;

	for (unsigned int i =0; i < kmer.length(); i++) {

		char c = kmer[i];
		char_map[c]++;

	}

	float entropy = 0;

	char nucs[] = { 'G', 'A', 'T', 'C'};

	for (unsigned int i = 0; i < 4; i++) {

		char nuc = nucs[i];
		int count = char_map[nuc];
		float prob = (float)count / kmer.length();

		if (prob > 0) {

			float val = prob * log(1/prob)/log(2.0f);
			entropy += val;

		}

	}

	return (entropy);
}


float compute_entropy (kmer_int_type kmer, unsigned int kmer_length){

	char count[] = {0, 0, 0, 0};

	for (unsigned int i = 0; i < kmer_length; i++) {

		int c = kmer & 3;
		kmer = kmer >> 2;
		count [c]++;

	}

	float entropy = 0;

	for (unsigned int i = 0; i < 4; i++) {

		float prob = (float)count[i] / kmer_length;

		if (prob > 0) {

			float val = prob * log(1/prob)/log(2.0f);
			entropy += val;

		}
	}

	return (entropy);

}


string revcomp (const string& kmer) {

	string revstring;

	for(int i = kmer.size() -1; i >= 0;i--) {

		char c = kmer[i];
		char revchar;

		switch (c) {

		case 'g':
			revchar = 'c';
			break;

		case 'G':
			revchar = 'C';
			break;

		case 'a':
			revchar = 't';
			break;

		case 'A':
			revchar = 'T';
			break;

		case 't':
			revchar = 'a';
			break;

		case 'T':
			revchar = 'A';
			break;

		case 'c':
			revchar = 'g';
			break;

		case 'C':
			revchar = 'G';
			break;

		default:
			revchar = 'N';
		}

		revstring += revchar;

	}

	return (revstring);

}


kmer_int_type revcomp_val (kmer_int_type kmer, unsigned int kmer_length) {

	kmer_int_type rev_kmer = 0;
	kmer = ~kmer;

	for(unsigned int i=0; i < kmer_length; i++) {

		int base = kmer & 3ll;
		rev_kmer = rev_kmer << 2;
		rev_kmer += base;
		kmer = kmer >> 2;

	}

	return rev_kmer;

}


kmer_int_type get_DS_kmer_val (kmer_int_type kmer, unsigned int kmer_length) {

	kmer_int_type rev_kmer = revcomp_val (kmer, kmer_length);

	if(rev_kmer > kmer)
		kmer = rev_kmer;

	return (kmer);
}

