#include "ReadUtility.h"
#include "GeneralSet.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>


using namespace std;

const int MAX_STR = 1024;


void load_reads(string file, bool rev) {

	
	fstream in;
	in.open(file.c_str(), fstream::in);

	if (!in.is_open()) {

		cout << "Error! Can't open file " << file << endl;
		exit(1);

	}

	cout << "Loading reads from file " << file << " ..." << endl;

	char temp[MAX_STR];

	while (!in.eof()) {

		in.getline(temp, MAX_STR);

		if(temp[0] == '>') {

			in.getline(temp, MAX_STR);
			string read(temp);			
			if (read.length() < g_read_length){
				for (int i = 0; i < g_read_length - read.length(); i++)
					read = read + "N";
			}
			if (read.length() > g_read_length){
				read = read.substr(0,g_read_length);
			}			
            if (rev)
                read = revcomp(read);
			data.push_back(read);			
			data_tag.push_back(-1);	
			
		}

	}

	cout <<"Total load " << data.size() << " reads!" << endl;
	in.close();

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
