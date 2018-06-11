#include "ReadUtility.h"
#include "GeneralSet.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>


using namespace std;

const int MAX_STR = 1024;


void load_reads(string file, vector<string>& input_data, bool rev) {
	
	fstream in;
	in.open(file.c_str(), fstream::in);

	if (!in.is_open()) {

		cout << "Error! Can't open file " << file << endl;
		exit(1);

	}

	cout << "Loading reads from file " << file << " ..." << endl;
	const int read_len = g_read_length;
	char temp[read_len];
	string read;
	read.reserve(read_len);
	in.getline(temp, MAX_STR);
	in.getline(temp, MAX_STR);
	read = temp;
	while (!in.eof()) {

		in.getline(temp, MAX_STR);	

		if(temp[0] == '>') {
		
            		if (rev)
               			read = revcomp(read);
			string (read).swap(read);
			input_data.push_back(read);			
			in.getline(temp, MAX_STR);
			read = temp;			
		} else {
			read = read + temp;		
		}

	}

	if (rev){
//cout << read << endl;
		read = revcomp(read);
//cout << read << endl;
	}
	string (read).swap(read);
	input_data.push_back(read);
	in.close();

}

void delete_error_reads(vector<string>& input_data) {

	cout << "Delete error reads..." << endl;
	time_t s_time = time(NULL);

	vector<string>(input_data).swap(data);
	data.clear();
	const size_t data_size = input_data.size();
	vector<bool> xxx(data_size, true);
	int max_read_id = input_data.size() / 2;
	for (int i = 0; i < input_data.size(); ++i) {
		const string& read = input_data[i];
		if (contains_non_gatc(read) || read.length() < g_read_length) {
			xxx[i] = false;
			int mate_id;
			if (i >= max_read_id)
				mate_id = i - max_read_id;
			else
				mate_id = i + max_read_id;
			xxx[mate_id] = false;
		} else {
			if (read.length() > g_read_length)
				input_data[i] = read.substr(0, g_read_length);
		}
	}
	for (int i = 0; i < input_data.size(); ++i) {
		if (xxx[i])
			data.push_back(input_data[i]);
	}

	vector<string>().swap(input_data);
	vector<string>(data).swap(data);
	const int data_size1 = data.size();
	vector<int> xxx1(data_size1, -5);
	vector<int>(xxx1).swap(data_tag);
	vector<int>().swap(xxx1);
	vector<bool>().swap(xxx);
	time_t e_time = time(NULL);
	cout << "Now, there are total of " << data.size() << " reads (elapsed time: " << (e_time - s_time) << " s)" << endl;
}



int count_reads_file(string file){
	
	
	fstream in;
	in.open(file.c_str(), fstream::in);

	if (!in.is_open()) {

		cout << "Error! Can't open file " << file << endl;
		exit(1);

	}

	cout << "Count reads from file " << file << " ..." << endl;

	char temp[MAX_STR];
	int total_reads = 0;

	while (!in.eof()) {
		in.getline(temp, MAX_STR);
		if(temp[0] == '>') {

			++total_reads;
		}

	}

	in.close();
	cout << "File " << file << " contains " << total_reads << "reads." << endl; 
	return total_reads;


}



read_int_type get_read_int_type(const string& read) {

	const int vec_len = g_read_length / 32;
	read_int_type read_int(vec_len+1);
	int i = 0;
	while ( i < vec_len) {
		const string& kmer = read.substr(i*32, 32);
		read_int[i] = kmer_to_int(kmer);		
		++i;
	}
	int remain_len = read.length() - vec_len*32;
	if (remain_len > 0) {
		const string& kmer = read.substr(read.length()-remain_len);
		read_int[i] = kmer_to_int(kmer);
	}
	
	return read_int;

}


void get_read_int_type(const string& read, read_int_type& read_int, const int vec_len) {

	int i = 0;
	while ( i < vec_len) {
		const string& kmer = read.substr(i*32, 32);
		read_int[i] = kmer_to_int(kmer);		
		++i;
	}
	int remain_len = read.length() - vec_len*32;
	if (remain_len > 0) {
		const string& kmer = read.substr(read.length()-remain_len);
		read_int[i] = kmer_to_int(kmer);
	}

}


/*

//functions in class InfoReadset
void InfoReadset::add(size_t i) {

	count++;
	readset.push_back(i);   //change by Jane

}


vector<size_t> InfoReadset::get_readset() {

	return readset;

}


size_t InfoReadset::get_count() {

	return count;

}


void InfoReadset::clear() {

	count = 0;
	readset.clear();

}
*/
