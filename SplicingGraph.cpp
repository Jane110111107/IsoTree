// This file was modified from a file named splicing_graph.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpacker Software LICENSE.
#include "SplicingGraph.h"
#include<iostream>

using namespace std;



bool is_seq_similar(const string& str1, const string& str2, char mode, double sim_error) {

	int mismatch = 0;
	int length = str1.length() > str2.length() ? str2.length(): str1.length();

	if (length == 0) {

		if (str1.empty() && str2.empty()) 
			return true;
		else 
			return false;
		
	}

	if (mode == 'F') {

		for (int i = 0; i < length; i++) {

			if (str1[i] != str2[i])
				mismatch++;

		}

	} else {

		for (int i = 0; i < length; i++) {

			if (str1[str1.length()-i-1] != str2[str2.length()-i-1]) 
				mismatch++;

		}

	}

	if ((float)mismatch/length < sim_error) 
		return true;
	else
		return false;

}


bool is_aligned(const string& str1, const string& str2, char mode, int tag) {

	int mismatch = 0;
	int length = str1.length() > str2.length()? str2.length(): str1.length();

	if (mode == 'F') {	
		for(int i = 0; i < length; i++) {	
			if(str1[i] != str2[i])						
				mismatch++;
		}
	} else {
		for (int i = 0; i < length; i++) {
			if (str1[str1.length()-i-1] != str2[str2.length()-i-1])
				mismatch++;
		}
	}
	return (mismatch <= tag);

}


bool compatible(const string& str1, const string& str2) {

	for(unsigned int i = 0; i <= str2.length()-g_kmer_length; i++) {

		const string& kmer = str2.substr(i,g_kmer_length);
		string::size_type start = str1.find(kmer);

		if(start != string::npos) {

			if (start > i)
				return is_aligned(str1.substr(start-i),str2,'F',2);
			else
				return is_aligned(str1,str2.substr(i-start),'F',2);

		}

	}
	
	return false;

}


SplicingGraph::SplicingGraph() {

	node_sum = 0;

}


size_t SplicingGraph::get_node_sum() {

	return node_sum;

}


void  SplicingGraph::set_parents() {

	for (int i = 0; i < node_sum; i++) {

		if (!node_set[i].parents.empty())
			node_set[i].parents.clear();

	}

	for (size_t i = 0; i < node_sum; i++) {

		vector<node_id>::iterator it;

		for (it = node_set[i].children.begin(); it != node_set[i].children.end(); it++)
			node_set[*it].add_parent(i);

	}

}



int SplicingGraph::add_node(Node& node) {

	node_set.push_back(node);
	return (node_sum++);

}



bool SplicingGraph::build(KmerHash& kmer_hash, int seed) {

	data_tag[seed] = -3;
	string seed_str = data[seed];
	list<int> map_reads;	
	cout << " Get trunk..." << endl;	
	const string& left = reverse_extend(kmer_hash, seed_str, map_reads);
	map_reads.push_back(seed);
	const string& right = forward_extend(kmer_hash, seed_str, map_reads);
	string trunk = left + right.substr(g_read_length);


	if(static_cast<int>(trunk.length()) < g_min_trunk_length) {
		cout << "The trunk is too short!" << endl;
		list<int>::iterator its;
		for (its = map_reads.begin(); its != map_reads.end(); its++) 
			data_tag[*its] = -1;		
		return false;
	}
	

	Node node(trunk);
	int p = add_node(node);	
	refine_forward(kmer_hash, p);
	refine_reverse(kmer_hash, p);

	list<int>::iterator its;
	for (its = map_reads.begin(); its != map_reads.end(); its++) 
		data_tag[*its] = -1;

	set_reads_pos_in_node(kmer_hash, p);
///*
	int node_len = node_set[p].sequence.length();
	vector<float> node_cov(node_len, 0.0);
	for (int i = 0; i < node_set[p].cov_reads.size(); i++)
		for (int j = node_set[p].cov_reads[i].second; j < node_set[p].cov_reads[i].second+g_read_length; j++)
			node_cov[j] = node_cov[j] + 1.0;

	bool first_reverse = true;
	bool first_forward = true;
	int reverse_start = 0;
	int reverse_end = 0;
	int forward_start = 0;
	int forward_end = 0;
	bool forward_tag = false;
	bool reverse_tag = false;
	vector<pair<int, int> > candi_reverse;
	vector<pair<int, int> > candi_forward;
	if (node_cov.size() > 2*g_read_length){
	for (int i = g_read_length; i < node_cov.size()-2*g_read_length; i++) {

		if (node_cov[i] / node_cov[i+g_min_same_len-1] < 0.5 || node_cov[i+g_min_same_len-1] - node_cov[i] > 100){

			if (forward_tag) {
				pair<int, int> fp;
				fp.first = forward_start;
				fp.second = forward_end;
				candi_forward.push_back(fp);
				//node_set[p].forward_check_pos.push_back(fp);
				first_forward = true;
				forward_tag = false;
			}

			if (first_reverse) {
				reverse_start = i;
				first_reverse = false;
			}
			reverse_end = i + g_min_same_len-1;
			reverse_tag = true;
		
		} else {

			if (node_cov[i+g_min_same_len-1] / node_cov[i] < 0.5 || node_cov[i] - node_cov[i+g_min_same_len-1] > 100) {

				if (reverse_tag) {
					pair<int,int> rp;
					rp.first = reverse_start;
					rp.second = reverse_end;
					candi_reverse.push_back(rp);
					//node_set[p].reverse_check_pos.push_back(rp);
					first_reverse = true;
					reverse_tag = false;
				}

				if (first_forward) {
					forward_start = i;
					first_forward = false;
				}
				forward_end = i + g_min_same_len - 1;
				forward_tag = true;
			}else {

				if (forward_tag) {
					pair<int, int> fp;
					fp.first = forward_start;
					fp.second = forward_end;
					candi_forward.push_back(fp);
					//node_set[p].forward_check_pos.push_back(fp);
					first_forward = true;
					forward_tag = false;
				}

				if (reverse_tag) {
					pair<int,int> rp;
					rp.first = reverse_start;
					rp.second = reverse_end;
					candi_reverse.push_back(rp);
					//node_set[p].reverse_check_pos.push_back(rp);
					first_reverse = true;
					reverse_tag = false;
				}
				
			}

		}

	}
	if (candi_forward.size() > 0) {
		int n = 1;
		pair <int, int> current_pos = candi_forward[0];
		
		while (n < candi_forward.size()) {
			if (candi_forward[n].first <= current_pos.second) {
				current_pos.second = candi_forward[n].second;
			} else {
				node_set[p].forward_check_pos.push_back(current_pos);
				current_pos = candi_forward[n];
			}		
			n = n + 1;
		}

		node_set[p].forward_check_pos.push_back(current_pos);
	}

	if (candi_reverse.size() > 0) {
		int n = 1;
		pair<int, int> current_pos = candi_reverse[0];
		while(n < candi_reverse.size()) {
			if (candi_reverse[n].first < current_pos.second)
				current_pos.second = candi_reverse[n].second;
			else {
				node_set[p].reverse_check_pos.push_back(current_pos);
				current_pos = candi_reverse[n];
			}
			n = n + 1;
		}
		node_set[p].reverse_check_pos.push_back(current_pos);
	}

	}
//*/
//set_check_pos(p);
	//
	if (g_mode == 1) {
		cout << "branch check and extend..." << endl;
		//branch_check_and_extend(kmer_hash, p);
		branch_extend_by_coverage(kmer_hash);
		cout << "refine graph..." << endl;
		refine_graph(kmer_hash);
		//trim_graph(kmer_hash);
	}

	return true;
/*
	if ((node_sum == 1 && (int)node_set[0].sequence.length() > g_min_transcript_length) 
		|| (node_sum > 1 && (int)get_total_amount_of_kmers() > g_min_kmers_per_graph)) 				
			return true;
	else {
		return false;
	}
*/
}



void SplicingGraph::refine_graph(KmerHash& kmer_hash) {

	for (int i = 0; i < node_sum; i++) {

		if (node_set[i].children.size() == 0)
			refine_forward(kmer_hash, i);
		
		if (node_set[i].parents.size() == 0)
			refine_reverse(kmer_hash, i);

		if (node_set[i].parents.size() == 0 || node_set[i].children.size() == 0)
			set_reads_pos_in_node(kmer_hash, i);

	}

}




void SplicingGraph::refine_forward(KmerHash& kmer_hash, node_id p) {

	if (p < 0 || p >= node_sum)
		return;
	if (node_set[p].sequence.length() <= g_read_length)
		return;
cout << "refine forward..." << endl;
	//set_reads_pos_in_node(kmer_hash, p);
	bool is_add = false;
	int check_len = g_refine_check_len;
	if (node_set[p].sequence.length() < g_refine_check_len + g_read_length)
		check_len = node_set[p].sequence.length()-g_read_length;

	for (int i = node_set[p].sequence.length()-g_read_length-1; i > node_set[p].sequence.length()-g_read_length-check_len; i--) {//forward check
		string seed_read = node_set[p].sequence.substr(i, g_read_length);
		list<int> map_reads;
		string forward_seq = forward_extend(kmer_hash, seed_read, map_reads);
		if (forward_seq.length() > g_min_exon_length + g_read_length && forward_seq.length() > node_set[p].sequence.length()-i && pair_support(map_reads, p)) {	
			string recover_str;
			if (static_cast<int>(i-g_read_length) > 0)
				recover_str = node_set[p].sequence.substr(i-g_read_length);//
			else
				recover_str = node_set[p].sequence;
			set_reads_tag(kmer_hash, recover_str, -1);//		
			node_set[p].sequence = node_set[p].sequence.substr(0, i) + forward_seq;	
			//set_reads_pos_in_node(kmer_hash, p);
			is_add = true;		
			break;
			
		} else {
			list<int>::iterator it;
			for(it = map_reads.begin(); it != map_reads.end(); it++)
				data_tag[*it] = -1;
		
		}
	}

	if (is_add)
		refine_forward(kmer_hash, p);

}



void SplicingGraph::refine_reverse(KmerHash& kmer_hash, node_id p) {

	if (p < 0 || p >= node_sum)
		return;
	if (node_set[p].sequence.length() <= g_read_length)
		return;
cout << "refine reverse..." << endl;
	//set_reads_pos_in_node(kmer_hash, p);
	bool is_add = false;
	int check_len = g_refine_check_len;
	if (node_set[p].sequence.length() < g_refine_check_len + g_read_length)
		check_len = node_set[p].sequence.length()-g_read_length;

	for (int i = 1; i < check_len; i++){//reverse check
		string seed_read = node_set[p].sequence.substr(i, g_read_length);		
		list<int> map_reads;
		string reverse_seq = reverse_extend(kmer_hash, seed_read, map_reads);
		if (reverse_seq.length() > g_min_exon_length + g_read_length && reverse_seq.length() > i + g_read_length && pair_support(map_reads, p)) {			
			string recover_str;
			if (i+2*g_read_length-1 < node_set[p].sequence.length())
				recover_str = node_set[p].sequence.substr(0, i+2*g_read_length-1);
			else
				recover_str = node_set[p].sequence;
			set_reads_tag(kmer_hash, recover_str, -1);//
			node_set[p].sequence = reverse_seq + node_set[p].sequence.substr(i+g_read_length);
			//set_reads_pos_in_node(kmer_hash, p);
			is_add = true;			
			break;

		}else {
			list<int>::iterator it;
			for(it = map_reads.begin(); it != map_reads.end(); it++)
				data_tag[*it] = -1;
		}

	}
	
	if (is_add)
		refine_reverse(kmer_hash, p);

}


bool SplicingGraph::pair_support(list<int> map_reads, node_id p) {

	if (p < 0 || p >= node_sum)
		return false;

	if (!g_is_paired_end)
		return true;

	int support = 0;
	list<int>::iterator it;
	for (it = map_reads.begin(); it != map_reads.end(); it++) {
		int mate_id;
		if (*it >= max_read_id)
			mate_id = *it - max_read_id;
		else
			mate_id = *it + max_read_id;
		if (data_tag[mate_id] == p) {
			support++;
			if (support >= g_min_reads_support_branch)
				break;
		}
	}

	if (support >= g_min_reads_support_branch)
		return true;
	else
		return false;

}


void SplicingGraph::get_reads(KmerHash& kmer_hash, string& sequence, list<int>& map_reads, char tag, int color_tag) {

	if (sequence.length() < g_read_length)
		return;

	if (tag == 'F') { //map_reads.push_front(*a)

		for (int i = sequence.length() - g_read_length; i >= 0 ; i--) {
			string real_read = sequence.substr(i,g_read_length);
			kmer_int_type a_kmer = kmer_to_int(sequence.substr(i, g_kmer_length)); 
			kmer_int_type b_kmer = kmer_to_int(sequence.substr(i+g_read_length-g_kmer_length, g_kmer_length));
			vector<int> a_readset = kmer_hash[a_kmer];
			vector<int> b_readset = kmer_hash[b_kmer];
			vector<int>::iterator a = a_readset.begin();
			vector<int>::iterator b = b_readset.begin();
			while (a != a_readset.end() && b != b_readset.end()) {

				if (*a == *b) {
					if (data_tag[*a] != -2) {
						string read = data[*a];				
						if (read == real_read)	{	
							data_tag[*a] = color_tag;
							map_reads.push_front(*a);
						}
						else {
							string read_rev = revcomp(read);
							if (g_is_paired_end && g_double_stranded_mode && read_rev == real_read){
								data_tag[*a] = color_tag;
								map_reads.push_front(*a);
							}
						}
					}
					a++;
					b++;

				} else {

					if (*a > *b)
						b++;
					else
						a++;
	
				}

			}
		}

	} else { //map_reads.push_fback(*a)

		for (int i = 0; i <= sequence.length() - g_read_length; i++) {

			string real_read = sequence.substr(i,g_read_length);
			kmer_int_type a_kmer = kmer_to_int(sequence.substr(i, g_kmer_length)); 
			kmer_int_type b_kmer = kmer_to_int(sequence.substr(i+g_read_length-g_kmer_length, g_kmer_length));
			vector<int> a_readset = kmer_hash[a_kmer];
			vector<int> b_readset = kmer_hash[b_kmer];
			vector<int>::iterator a = a_readset.begin();
			vector<int>::iterator b = b_readset.begin();
			while (a != a_readset.end() && b != b_readset.end()) {

				if (*a == *b) {
					if (data_tag[*a] != -2) {
						string read = data[*a];				
						if (read == real_read) {				
							map_reads.push_back(*a);
							data_tag[*a] = color_tag;
						} else {
							string read_rev = revcomp(read);
							if (g_is_paired_end && g_double_stranded_mode && read_rev == real_read){
								map_reads.push_back(*a);
								data_tag[*a] = color_tag;
							}
						}
					}
					a++;
					b++;

				} else {

					if (*a > *b)
						b++;
					else
						a++;
	
				}

			}


		}


	}

}



string SplicingGraph::forward_extend(KmerHash& kmer_hash, string seed_contig, list<int>& map_reads, pair<int, int>& reach_id, int color_tag, bool reUse) {

	//cout << "forward extend..." << endl;
	string contig = seed_contig;
	bool stop_tag = false;
	if (contig.length() < g_read_length)
		return contig;
	map<string, vector<int> > candidates;
	map<string, int> candi_pairs;
	while (!stop_tag) {

		bool has_reach = false;
		reach_id.first = -1;
		reach_id.second = -1;
		string forward_kmer = contig.substr(contig.length()-g_kmer_length, g_kmer_length);
		kmer_int_type kmer_int = kmer_to_int(forward_kmer);
		vector<int> seed_readset = kmer_hash[kmer_int];
		//
		if (seed_readset.size() == 0) {
			stop_tag = true;
			continue;
		}
				
		for (int i = g_max_same_len; i >= g_min_same_len; i--) {
			string first_kmer = contig.substr(contig.length()-i,g_kmer_length);
			kmer_int_type tag_kmer = kmer_to_int(first_kmer);	 		
	        	vector<int> tag_readset = kmer_hash[tag_kmer];
			string ref_seq = contig.substr(contig.length()-i);
	       		vector<int>::iterator a = seed_readset.begin();
			vector<int>::iterator b = tag_readset.begin();
			candidates.clear();
			candi_pairs.clear();
			bool stop_for = false;
	        while (a != seed_readset.end() && b != tag_readset.end()) {
		  
				if (*a == *b) {

					if ((!has_reach) && data_tag[*a] >= 0 && contig.length() > seed_contig.length()) {
						reach_id.first = *a;
						reach_id.second = i;
						has_reach = true;
					}
					int pair_support = 0;
					if (*a >= max_read_id && (data_tag[*a-max_read_id] >= 0 || data_tag[*a-max_read_id] == color_tag))
						pair_support = 1;
					if (*a < max_read_id && (data_tag[*a+max_read_id] >= 0 || data_tag[*a+max_read_id] == color_tag))
						pair_support = 1;

					if ((!reUse) && (data_tag[*a] == -1)) {

						if (is_aligned(ref_seq, data[*a], 'F')){ 
							string candi = data[*a].substr(i);
							candidates[candi].push_back(*a);
							if (candidates[candi].size() == 1) {
								candi_pairs[candi] = pair_support;
							} else {
								if (candi_pairs[candi] == 0)
									candi_pairs[candi] = pair_support;
							}

						} else if (g_is_paired_end && g_double_stranded_mode) {

							string rev_read = revcomp(data[*a]);
							if (is_aligned(ref_seq, rev_read, 'F')) {
								string candi = rev_read.substr(i);
								candidates[candi].push_back(*a);
								if (candidates[candi].size() == 1) {
									candi_pairs[candi] = pair_support;
								} else {
									if (candi_pairs[candi] == 0)
										candi_pairs[candi] = pair_support;
								}
							}
						}
						
					} else if (reUse) {

						if (is_aligned(ref_seq, data[*a], 'F')){ 
							string candi = data[*a].substr(i);
							candidates[candi].push_back(*a);
							if (candidates[candi].size() == 1) {
								candi_pairs[candi] = pair_support;
							} else {
								if (candi_pairs[candi] == 0)
									candi_pairs[candi] = pair_support;
							}

						} else if (g_is_paired_end && g_double_stranded_mode) {

							string rev_read = revcomp(data[*a]);
							if (is_aligned(ref_seq, rev_read, 'F')) {
								string candi = rev_read.substr(i);
								candidates[candi].push_back(*a);
								if (candidates[candi].size() == 1) {
									candi_pairs[candi] = pair_support;
								} else {
									if (candi_pairs[candi] == 0)
										candi_pairs[candi] = pair_support;
								}
							}
						}
					}
								
					a++;	
					b++;
	
				} else {

					if (*a > *b)		
						b++;		
					else 			
						a++;	
				}
			}
			
			if (stop_for)
				break;
			if (has_reach && candidates.size() == 0) {
				stop_tag = true;
				break;
			}
			if (candidates.size() > 0) {
				map<string, vector<int> >::iterator it = candidates.begin();
				int max_candi = it->second.size();
				string extern_str = it -> first;
				string extern_str1;
				
				for (; it !=candidates.end(); it++) {
					if (it -> second.size() > max_candi) {
						max_candi = it -> second.size();
						extern_str = it -> first;
						if (candi_pairs[extern_str] == 1)
							extern_str1 = extern_str;
					}
				}
				

				if (extern_str1.length() > 0){
					contig = contig + extern_str1;
					for (int j = 0; j < candidates[extern_str1].size(); j++) {
						data_tag[candidates[extern_str1][j]] = color_tag;
						map_reads.push_back(candidates[extern_str1][j]);
					}
				} else{
					/*
					srand((unsigned)time(NULL));
					int rand_pos = rand() % candidates.size();
					it = candidates.begin();
					for(int k = 0; k < rand_pos; k++)
						it++;
					extern_str = it->first;
					*/
					contig = contig + extern_str;
					for (int j = 0; j < candidates[extern_str].size(); j++) {
						data_tag[candidates[extern_str][j]] = color_tag;
						map_reads.push_back(candidates[extern_str][j]);
					}
				}
				break;
			} else if (i == g_min_same_len)
				stop_tag = true;
		}

	}

	return contig;

}


string SplicingGraph::forward_extend(KmerHash& kmer_hash, string seed_contig, list<int>& map_reads, int color_tag, bool reUse) {

	//cout << "forward extend..." << endl;
	string contig = seed_contig;
	bool stop_tag = false;
	if (contig.length() < g_read_length)
		return contig;
	map<string, vector<int> > candidates;
	map<string, int> candi_pairs;
	while (!stop_tag) {


		string forward_kmer = contig.substr(contig.length()-g_kmer_length, g_kmer_length);
		kmer_int_type kmer_int = kmer_to_int(forward_kmer);
		vector<int> seed_readset = kmer_hash[kmer_int];

		if (seed_readset.size() == 0) {
			stop_tag = true;
			continue;
		}
		//
/*		
		if (g_is_paired_end){
//cout <<"if (g_is_paired_end){"<<endl;
			for (int i = 0; i < seed_readset.size(); i++) {
				if (data_tag[seed_readset[i]] != -1)
					continue;
				int mate_id;
				if (seed_readset[i] >= max_read_id)
					mate_id = seed_readset[i] - max_read_id;
				else
					mate_id = seed_readset[i] + max_read_id;
				if (data_tag[mate_id] == color_tag || data_tag[mate_id] >= 0) {
					string read = data[seed_readset[i]];
					string::size_type start = read.find(forward_kmer);
					if (start == string::npos && g_double_stranded_mode) {
						read = revcomp(read);
						start = read.find(forward_kmer);
					}
					if (start == string::npos|| start == g_read_length-g_kmer_length)
						continue;
					//if (is_aligned(read.substr(0, start+g_kmer_length), contig.substr(contig.length()-start-g_kmer_length))) {
					if (read.substr(0, start+g_kmer_length) == contig.substr(contig.length()-start-g_kmer_length)) {
//cout << "if =="<<endl;
						string add_reads_seq = contig.substr(contig.length()-g_read_length + 1) + read.substr(start+g_kmer_length);
						//data_tag[seed_readset[i]] = color_tag;
						get_reads(kmer_hash, add_reads_seq, map_reads, 'B', color_tag);
						contig = contig + read.substr(start+g_kmer_length);												
						continue_tag = true;
						break;
					}
				}

			}
		}

		if (continue_tag)
			continue;
//cout << "if-continue"<<endl;
*/
		

		for (int i = g_max_same_len; i >= g_min_same_len; i--) {
			string first_kmer = contig.substr(contig.length()-i,g_kmer_length);
			kmer_int_type tag_kmer = kmer_to_int(first_kmer);	 
	        vector<int> tag_readset = kmer_hash[tag_kmer];
			string ref_seq = contig.substr(contig.length()-i);
	        vector<int>::iterator a = seed_readset.begin();
			vector<int>::iterator b = tag_readset.begin();
			candidates.clear();	
			candi_pairs.clear();
	        while (a != seed_readset.end() && b != tag_readset.end()) {
		  
				if (*a == *b) {

					int pair_support = 0;
					if (*a >= max_read_id && (data_tag[*a-max_read_id] >= 0 || data_tag[*a-max_read_id] == color_tag))
						pair_support = 1;
					if (*a < max_read_id && (data_tag[*a+max_read_id] >= 0 || data_tag[*a+max_read_id] == color_tag))
						pair_support = 1;


					if ((!reUse) && (data_tag[*a] == -1)) {

						if (is_aligned(ref_seq, data[*a], 'F')){ 
							string candi = data[*a].substr(i);
							candidates[candi].push_back(*a);
							if (candidates[candi].size() == 1) {
								candi_pairs[candi] = pair_support;
							} else {
								if (candi_pairs[candi] == 0)
									candi_pairs[candi] = pair_support;
							}

						} else if (g_is_paired_end && g_double_stranded_mode) {

							string rev_read = revcomp(data[*a]);
							if (is_aligned(ref_seq, rev_read, 'F')) {
								string candi = rev_read.substr(i);
								candidates[candi].push_back(*a);
								if (candidates[candi].size() == 1) {
									candi_pairs[candi] = pair_support;
								} else {
									if (candi_pairs[candi] == 0)
										candi_pairs[candi] = pair_support;
								}
							}
						}
						
					} else if (reUse) {

						if (is_aligned(ref_seq, data[*a], 'F')){ 
							string candi = data[*a].substr(i);
							candidates[candi].push_back(*a);
							if (candidates[candi].size() == 1) {
								candi_pairs[candi] = pair_support;
							} else {
								if (candi_pairs[candi] == 0)
									candi_pairs[candi] = pair_support;
							}

						} else if (g_is_paired_end && g_double_stranded_mode) {

							string rev_read = revcomp(data[*a]);
							if (is_aligned(ref_seq, rev_read, 'F')) {
								string candi = rev_read.substr(i);
								candidates[candi].push_back(*a);
								if (candidates[candi].size() == 1) {
									candi_pairs[candi] = pair_support;
								} else {
									if (candi_pairs[candi] == 0)
										candi_pairs[candi] = pair_support;
								}
							}
						}
					}
								
					a++;	
					b++;
	
				} else {

					if (*a > *b)		
						b++;		
					else 			
						a++;	
				}
			}

			if (candidates.size() > 0) {
				map<string, vector<int> >::iterator it = candidates.begin();
				int max_candi = it->second.size();
				string extern_str = it -> first;
				string extern_str1;
				for (; it !=candidates.end(); it++) {
					if (it -> second.size() > max_candi) {
						max_candi = it -> second.size();
						extern_str = it -> first;
						if (candi_pairs[extern_str] == 1)
							extern_str1 = extern_str;
					}
				}
				if (extern_str1.length() > 0) {
					contig = contig + extern_str1;
					for (int j = 0; j < candidates[extern_str1].size(); j++) {
						data_tag[candidates[extern_str1][j]] = color_tag;
						map_reads.push_back(candidates[extern_str1][j]);
					}
				}else{
					/*
					srand((unsigned)time(NULL));
					int rand_pos = rand() % candidates.size();
					it = candidates.begin();
					for(int k = 0; k < rand_pos; k++)
						it++;
					extern_str = it->first;
					*/
					contig = contig + extern_str;
					for (int j = 0; j < candidates[extern_str].size(); j++) {
						data_tag[candidates[extern_str][j]] = color_tag;
						map_reads.push_back(candidates[extern_str][j]);
					}
				}
				break;
			} else if (i == g_min_same_len)
				stop_tag = true;
		}

	}

	return contig;

}




string SplicingGraph::reverse_extend(KmerHash& kmer_hash,  string seed_contig, list<int>& map_reads,  pair<int, int>& reach_id, int color_tag, bool reUse) {

	//cout << "reverse extend..." << endl;
	string contig = seed_contig;
	if (contig.length() < g_read_length)
		return contig;

	bool stop_tag = false;
	map<string, vector<int> > candidates;
	map<string, int> candi_pairs;

	while (!stop_tag) {
		
		bool reach_tag = false;
		reach_id.first = -1;
		reach_id.second = -1;
		string reverse_kmer = contig.substr(0,g_kmer_length);
		kmer_int_type kmer_int = kmer_to_int(reverse_kmer);
		vector<int> seed_readset = kmer_hash[kmer_int];
		//
		if (seed_readset.size() == 0) {
			stop_tag = true;
			continue;
		}
		//
/*		
		if (g_is_paired_end){
//cout <<"if (g_is_paired_end){"<<endl;			
			for (int i = 0; i < seed_readset.size(); i++) {
				if (data_tag[seed_readset[i]] != -1)
					continue;
				int mate_id;
				if (seed_readset[i] >= max_read_id)
					mate_id = seed_readset[i] - max_read_id;
				else
					mate_id = seed_readset[i] + max_read_id;
				if (data_tag[mate_id] == color_tag || data_tag[mate_id] >= 0) {
					string read = data[seed_readset[i]];
					string::size_type start = read.find(reverse_kmer);
					if (start == string::npos && g_double_stranded_mode) {
						read = revcomp(read);
						start = read.find(reverse_kmer);
					}
					if (start == string::npos || start == 0)
						continue;
					//if (is_aligned(read.substr(start), contig.substr(0, g_read_length-start))) {
					if (read.substr(start)==contig.substr(0, g_read_length-start)) {
//cout << "if =="<<endl;
						string add_reads_seq = read.substr(0, start) + contig.substr(0, g_read_length-1);
						get_reads(kmer_hash, add_reads_seq, map_reads, 'F', color_tag);
                        //data_tag[seed_readset[i]] = color_tag;
						contig = read.substr(0, start) + contig;												
						continue_tag = true;
						break;
					}
				}

			}
		}

		if (continue_tag)
			continue;
//cout << "if-continue"<<endl;
*/

		bool stop_for = false;
		for (int i = g_max_same_len-g_kmer_length; i >=g_min_same_len-g_kmer_length; i--) {
			string second_kmer = contig.substr(i,g_kmer_length);		
			kmer_int_type tag_kmer = kmer_to_int(second_kmer);
	        vector<int> tag_readset = kmer_hash[tag_kmer];
			string ref_seq = contig.substr(0,i+g_kmer_length);
		    vector<int>::iterator a = seed_readset.begin();
			vector<int>::iterator b = tag_readset.begin();
			candidates.clear();
			candi_pairs.clear();
	        while (a != seed_readset.end() && b != tag_readset.end()) {
		  
				if (*a == *b) {

					if ((!reach_tag) && data_tag[*a] >= 0 && contig.length() > seed_contig.length()) {
						reach_id.first = *a;
						reach_id.second = i + g_kmer_length;
						reach_tag = true;
					}

					int pair_support = 0;
					if (*a >= max_read_id && (data_tag[*a-max_read_id] >= 0 || data_tag[*a-max_read_id] == color_tag))
						pair_support = 1;
					if (*a < max_read_id && (data_tag[*a+max_read_id] >= 0 || data_tag[*a+max_read_id] == color_tag))
						pair_support = 1;

					if ((!reUse) && (data_tag[*a] == -1))	{

						if (is_aligned(ref_seq, data[*a], 'R')) {
							string candi = data[*a].substr(0,g_read_length-i-g_kmer_length);
							candidates[candi].push_back(*a);
							if (candidates[candi].size() == 1) {
								candi_pairs[candi] = pair_support;
							} else {
								if (candi_pairs[candi] == 0)
									candi_pairs[candi] = pair_support;
							}
						} else if (g_is_paired_end && g_double_stranded_mode) {
							string rev_read = revcomp(data[*a]);
							if (is_aligned(ref_seq, rev_read, 'R')) {
							    string candi = rev_read.substr(0,g_read_length-i-g_kmer_length);
							    candidates[candi].push_back(*a);
								if (candidates[candi].size() == 1) {
									candi_pairs[candi] = pair_support;
								} else {
									if (candi_pairs[candi] == 0)
										candi_pairs[candi] = pair_support;
								}
							}
						}

					} else if (reUse) {

						if (is_aligned(ref_seq, data[*a], 'R')) {
							string candi = data[*a].substr(0,g_read_length-i-g_kmer_length);
							candidates[candi].push_back(*a);	
							if (candidates[candi].size() == 1) {
								candi_pairs[candi] = pair_support;
							} else {
								if (candi_pairs[candi] == 0)
									candi_pairs[candi] = pair_support;
							}
						} else if (g_is_paired_end && g_double_stranded_mode) {
							string rev_read = revcomp(data[*a]);
							if (is_aligned(ref_seq, rev_read, 'R')) {
							    string candi = rev_read.substr(0,g_read_length-i-g_kmer_length);
							    candidates[candi].push_back(*a);
								if (candidates[candi].size() == 1) {
									candi_pairs[candi] = pair_support;
								} else {
									if (candi_pairs[candi] == 0)
										candi_pairs[candi] = pair_support;
								}
							}
						}
					}						
								
					a++;	
					b++;
	
				} else {

					if (*a > *b)		
						b++;		
					else 			
						a++;	
				}
			}

			if (stop_for)
				break;
			if (reach_tag && candidates.size() == 0) {
				stop_tag = true;
				break;
			}
			if (candidates.size() > 0) {

				map<string, vector<int> >::iterator it = candidates.begin();
				int max_candi = it->second.size();
				string extern_str = it -> first;
				string extern_str1;
				for (; it !=candidates.end(); it++) {
					if (it -> second.size() > max_candi) {
						max_candi = it -> second.size();
						extern_str = it -> first;
						if (candi_pairs[extern_str] == 1)
							extern_str1 = extern_str;
					}
				}

				if (extern_str1.length() > 0) {
					contig = extern_str1 + contig;		
					for (int j = 0; j < candidates[extern_str1].size(); j++) {
						data_tag[candidates[extern_str1][j]] = color_tag;
						map_reads.push_front(candidates[extern_str1][j]);
					}
				}else {
					/*
					srand((unsigned)time(NULL));
					int rand_pos = rand() % candidates.size();
					it = candidates.begin();
					for(int k = 0; k < rand_pos; k++)
						it++;
					extern_str = it->first;
					*/
					contig = extern_str + contig;		
					for (int j = 0; j < candidates[extern_str].size(); j++) {
						data_tag[candidates[extern_str][j]] = color_tag;
						map_reads.push_front(candidates[extern_str][j]);
					}
				}
				break;

			} else if (i == g_min_same_len-g_kmer_length)
				stop_tag = true;

		} //for (int i = g_max_same_len-g_kmer_length; i >=g_min_same_len-g_kmer_length; i--)
	} //while (!stop_tag)
	
	return contig;

}



string SplicingGraph::reverse_extend(KmerHash& kmer_hash,  string seed_contig, list<int>& map_reads, int color_tag, bool reUse) {

	//cout << "reverse extend..." << endl;
	string contig = seed_contig;
	if (contig.length() < g_read_length)
		return contig;

	bool stop_tag = false;
	map<string, vector<int> > candidates;
	map<string, int> candi_pairs;
	while (!stop_tag) {

		string reverse_kmer = contig.substr(0,g_kmer_length);
		kmer_int_type kmer_int = kmer_to_int(reverse_kmer);
		vector<int> seed_readset = kmer_hash[kmer_int];
		//
		if (seed_readset.size() == 0) {
			stop_tag = true;
			continue;
		}
		//
/*		
		if (g_is_paired_end){
//cout << "if (g_is_paired_end){" << endl;			
			for (int i = 0; i < seed_readset.size(); i++) {
				if (data_tag[seed_readset[i]] != -1)
					continue;
				int mate_id;
				if (seed_readset[i] >= max_read_id)
					mate_id = seed_readset[i] - max_read_id;
				else
					mate_id = seed_readset[i] + max_read_id;
				if (data_tag[mate_id] == color_tag || data_tag[mate_id] >= 0) {
					string read = data[seed_readset[i]];
					string::size_type start = read.find(reverse_kmer);
					if (start == string::npos && g_double_stranded_mode) {
						read = revcomp(read);
						start = read.find(reverse_kmer);
					}
					if (start == string::npos || start == 0)
						continue;
					//if (is_aligned(read.substr(start), contig.substr(0, g_read_length-start))) {
					if (read.substr(start)==contig.substr(0, g_read_length-start)) {
//cout << "if =="<<endl;
						string add_reads_seq = read.substr(0, start) + contig.substr(0, g_read_length-1);
						get_reads(kmer_hash, add_reads_seq, map_reads, 'F', color_tag);
					 	//data_tag[seed_readset[i]] = color_tag;
						contig = read.substr(0, start) + contig;												
						continue_tag = true;
						break;
					}
				}

			}
		}

		if (continue_tag)
			continue;
//cout << "if-continue"<<endl;
*/
		for (int i = g_max_same_len-g_kmer_length; i >=g_min_same_len-g_kmer_length; i--) {
			string second_kmer = contig.substr(i,g_kmer_length);		
			kmer_int_type tag_kmer = kmer_to_int(second_kmer);
			vector<int> seed_readset = kmer_hash[kmer_int];
	        vector<int> tag_readset = kmer_hash[tag_kmer];
			string ref_seq = contig.substr(0,i+g_kmer_length);
		    vector<int>::iterator a = seed_readset.begin();
			vector<int>::iterator b = tag_readset.begin();
			candidates.clear();
			candi_pairs.clear();
	        while (a != seed_readset.end() && b != tag_readset.end()) {
		  
				if (*a == *b) {

					int pair_support = 0;
					if (*a >= max_read_id && (data_tag[*a-max_read_id] >= 0 || data_tag[*a-max_read_id] == color_tag))
						pair_support = 1;
					if (*a < max_read_id && (data_tag[*a+max_read_id] >= 0 || data_tag[*a+max_read_id] == color_tag))
						pair_support = 1;

					if ((!reUse) && (data_tag[*a] == -1))	{

						if (is_aligned(ref_seq, data[*a], 'R')) {
							string candi = data[*a].substr(0,g_read_length-i-g_kmer_length);
							candidates[candi].push_back(*a);	
							if (candidates[candi].size() == 1) {
								candi_pairs[candi] = pair_support;
							} else {
								if (candi_pairs[candi] == 0)
									candi_pairs[candi] = pair_support;
							}
						} else if (g_is_paired_end && g_double_stranded_mode) {
							string rev_read = revcomp(data[*a]);
							if (is_aligned(ref_seq, rev_read, 'R')) {
							    string candi = rev_read.substr(0,g_read_length-i-g_kmer_length);
							    candidates[candi].push_back(*a);
								if (candidates[candi].size() == 1) {
									candi_pairs[candi] = pair_support;
								} else {
									if (candi_pairs[candi] == 0)
										candi_pairs[candi] = pair_support;
								}
							}
						}

					} else if (reUse) {

						if (is_aligned(ref_seq, data[*a], 'R')) {
							string candi = data[*a].substr(0,g_read_length-i-g_kmer_length);
							candidates[candi].push_back(*a);
							if (candidates[candi].size() == 1) {
								candi_pairs[candi] = pair_support;
							} else {
								if (candi_pairs[candi] == 0)
									candi_pairs[candi] = pair_support;
							}
						} else if (g_is_paired_end && g_double_stranded_mode) {
							string rev_read = revcomp(data[*a]);
							if (is_aligned(ref_seq, rev_read, 'R')) {
							    string candi = rev_read.substr(0,g_read_length-i-g_kmer_length);
							    candidates[candi].push_back(*a);
								if (candidates[candi].size() == 1) {
									candi_pairs[candi] = pair_support;
								} else {
									if (candi_pairs[candi] == 0)
										candi_pairs[candi] = pair_support;
								}
							}
						}
					}						
								
					a++;	
					b++;
	
				} else {

					if (*a > *b)		
						b++;		
					else 			
						a++;	
				}
			}


			if (candidates.size() > 0) {

				map<string, vector<int> >::iterator it = candidates.begin();
				int max_candi = it->second.size();
				string extern_str = it -> first;
				string extern_str1;
				for (; it !=candidates.end(); it++) {
					if (it -> second.size() > max_candi) {
						max_candi = it -> second.size();
						extern_str = it -> first;
						if (candi_pairs[extern_str] == 1)
							extern_str1 = extern_str;
					}
				}

				if (extern_str1.length() > 0) {
					contig = extern_str1 + contig;
					for (int j = 0; j < candidates[extern_str1].size(); j++) {
						data_tag[candidates[extern_str1][j]] = color_tag;
						map_reads.push_front(candidates[extern_str1][j]);
					}
				}else{
					/*
					srand((unsigned)time(NULL));
					int rand_pos = rand() % candidates.size();
					it = candidates.begin();
					for(int k = 0; k < rand_pos; k++)
						it++;
					extern_str = it->first;
					*/
					contig = extern_str + contig;
					for (int j = 0; j < candidates[extern_str].size(); j++) {
						data_tag[candidates[extern_str][j]] = color_tag;
						map_reads.push_front(candidates[extern_str][j]);
					}
				}
				break;

			} else if (i == g_min_same_len-g_kmer_length)
				stop_tag = true;

		} //for (int i = g_max_same_len-g_kmer_length; i >=g_min_same_len-g_kmer_length; i--)
	} //while (!stop_tag)
	
	return contig;

}


bool SplicingGraph::is_circle(node_id p, node_id q, set<node_id>& checked) {

	bool flag = false;
	if (node_set[q].children.empty() || checked.find(q) != checked.end())
		return false;
	else
		checked.insert(q);

	 vector<node_id>::iterator it;
	 for (it = node_set[q].children.begin(); it != node_set[q].children.end(); it++) {

		 if ((*it) == p) {
			 flag = true;
			 break;
		 } else {
			 flag = is_circle(p, *it, checked);
			 if (flag)
				 break;
		 }

	 }

	 return flag;

}


void SplicingGraph::trim_graph(KmerHash& kmer_hash) {

	set_coverage_of_nodes();
	bool has_trimed = false;
	map<pair_t, double> edge_coverages;
	vector<pair_t> edges;
	get_coverage_of_edges(kmer_hash, edge_coverages, edges);
	map<node_id, double> total_out_coverages;
	map<node_id, double> total_in_coverages;


	for (int i = 0; i < edges.size(); i++) {

		int source = edges[i].first;
		if (total_out_coverages.find(source) == total_out_coverages.end())
			total_out_coverages[source] = edge_coverages[edges[i]];
		else
			total_out_coverages[source] += edge_coverages[edges[i]];

		int target = edges[i].second;
		if (total_in_coverages.find(target) == total_in_coverages.end())
			total_in_coverages[target] = edge_coverages[edges[i]];
		else
			total_in_coverages[target] += edge_coverages[edges[i]];

	}

	for (int i = 0; i < edges.size(); i++) {

		int source = edges[i].first;
		int target = edges[i].second;

		//if ((node_set[source].parents.empty() && node_set[source].children.size() == 1) || 
		//	(node_set[target].children.empty() && node_set[target].parents.size() == 1)) {

				string check_edge = get_edge_sequence(source, target);
				if ((int)check_edge.length() < g_read_length)
					continue;
				double e_cov = edge_coverages[edges[i]];
				double flanking_node_cov = node_set[source].node_coverage > node_set[target].node_coverage ? node_set[source].node_coverage : node_set[target].node_coverage;
				
				if (e_cov < g_min_junction_coverage ||
					e_cov < g_min_ratio_welds * flanking_node_cov ||
					e_cov < g_min_ratio_branch * total_out_coverages[source] ||
					e_cov < g_min_ratio_branch * total_in_coverages[target] ||
					(total_in_coverages.find(source) != total_in_coverages.end() && e_cov < g_min_ratio_in_out * total_in_coverages[source]) ||
					(total_out_coverages.find(target) != total_out_coverages.end() && e_cov < g_min_ratio_in_out * total_out_coverages[target])) {

						has_trimed = true;
						set_reads_tag(kmer_hash, check_edge, -1);
						node_set[source].delete_child(target);
						node_set[target].delete_parent(source);
				}

		//}

	}
	set_parents();
}



void SplicingGraph::set_coverage_of_nodes() {

	for (int i = 0; i < node_sum; i++) {

		if (node_set[i].sequence.length() == 0) {
			node_set[i].node_coverage = 0;
			continue;
		}

		node_set[i].node_coverage = static_cast<float> (node_set[i].cov_reads.size()) * g_read_length / node_set[i].sequence.length();

	}

}


void SplicingGraph::get_coverage_of_edges(KmerHash& kmer_hash, map<pair_t, double>& edge_coverages, vector<pair_t>& edges) {

	edges.clear();
	edge_coverages.clear();

	for (int i = 0; i < node_sum; i++) {

		if (node_set[i].children.empty())
			continue;

		for (int j = 0; j < (int)node_set[i].children.size(); j++)
			edges.push_back(pair_t(i, node_set[i].children[j]));

	}

	for (int i = 0; i < (int)edges.size(); i++) {

		int source = edges[i].first;
		int target = edges[i].second;
		const string& edge = get_edge_sequence(source, target);
		edge_coverages[edges[i]] = compute_coverage(kmer_hash,edge);

	}

}


void SplicingGraph::get_coverage_of_edges(KmerHash& kmer_hash, map<pair_t, double>& edge_coverages) {

	vector<pair_t> edges;

	for (int i = 0; i < node_sum; i++) {

		if (node_set[i].children.empty())
			continue;

		for (int j = 0; j < (int)node_set[i].children.size(); j++)
			edges.push_back(pair_t(i, node_set[i].children[j]));

	}

	for (int i = 0; i < (int)edges.size(); i++) {

		int source = edges[i].first;
		int target = edges[i].second;
		const string& edge = get_edge_sequence(source, target);
		edge_coverages[edges[i]] = compute_coverage(kmer_hash, edge);

	}

}


float SplicingGraph::compute_coverage(KmerHash& kmer_hash, const string& sequence) {


	float total = 0.0;

	if (sequence.length() < g_read_length)
		return 0.0;
	
	if (contains_non_gatc(sequence))
		return 0.0;
	for (int i = 0; i <= sequence.length()-g_read_length; i++) {
		string real_read = sequence.substr(i,g_read_length);
		kmer_int_type a_kmer = kmer_to_int(sequence.substr(i, g_kmer_length)); 
		kmer_int_type b_kmer = kmer_to_int(sequence.substr(i+g_read_length-g_kmer_length, g_kmer_length));
		vector<int> a_readset = kmer_hash[a_kmer];
		vector<int> b_readset = kmer_hash[b_kmer];
		vector<int>::iterator a = a_readset.begin();
		vector<int>::iterator b = b_readset.begin();
		while (a != a_readset.end() && b != b_readset.end()) {

			if (*a == *b) {
				if (data_tag[*a] != -2) {
					string read = data[*a];
					string read_rev = revcomp(read);
					if (read == real_read)				
						total++;
					else {
				   		if (g_is_paired_end && g_double_stranded_mode && read_rev == real_read)
					    	total++;
					}
				}
				a++;
				b++;

			} else {

				if (*a > *b)
					b++;
				else
					a++;

			}

		}

	}

	total = total * g_read_length / sequence.length();
	return total;

}


string SplicingGraph::get_edge_sequence(node_id s, node_id t) {

	if (s < 0 || s >= node_sum)
		return "";
	if (t < 0 || t >= node_sum)
		return "";

	int length = g_read_length - 1;
	int start = (int)node_set[s].sequence.length() > length ? static_cast<int>(node_set[s].sequence.length())-length : 0;
	string edge = node_set[s].sequence.substr(start);

	while ((int)edge.length() < length) {

		if (node_set[s].parents.size() >= 1) {

			int parent = node_set[s].parents[0];
			int remain = length - edge.length();
			int start_p = (int)node_set[parent].sequence.length() > remain ? node_set[parent].sequence.length()-remain : 0;
			edge = node_set[parent].sequence.substr(start_p) + edge;
			s = parent;

		} else {

			break;

		}

	}
	int pre_len = edge.length();
	if ((int)node_set[t].sequence.length() < length)
		edge = edge + node_set[t].sequence;
	else
		edge = edge + node_set[t].sequence.substr(0, length);

	while ( (int)edge.length() < pre_len+length) {
		
		if (node_set[t].children.size() >= 1) {

			int child = node_set[t].children[0];
			int remain =pre_len + length - edge.length();
			if (node_set[child].sequence.length() >= remain)
				edge = edge + node_set[child].sequence.substr(0,remain);
			else 
				edge = edge + node_set[child].sequence;
			t = child;

		} else {

			break;

		}

	}

	// refine
	//int edge_len = edge.length();

	//if (edge_len > g_read_length + 8)
		//edge = edge.substr(2, edge_len-4);
	//else if (edge_len > g_read_length + 4)
		//edge = edge.substr(1, edge_len-2);

	return edge;

}


void SplicingGraph::get_transcripts(KmerHash& kmer_hash, vector<pair<string,float> >& transcripts) {
	cout << "get_transcripts..." << endl;
	set_coverage_of_nodes();
	TreeStruct isotree;
	TreeNode root_node(node_order[node_order.size()-1],0);
 	int root_id = isotree.add_node(root_node);	
	vector<pair<node_id,float> > leaf;
	vector<int> tree_leaf_id;
	vector<vector<int> > trans_node_vec;
	vector<int> dealed;

	for (int i = 0; i < node_set[node_order[node_order.size()-1]].children.size(); i++) {

		node_id root_child = node_set[node_order[node_order.size()-1]].children[i];

		if (node_set[root_child].node_coverage > g_min_ratio_non_error) {
		
			float node_in = 0;
			float node_out = 0;

			for (int j = 0; j < node_set[root_child].parents.size(); j++) {
				const string& edge_in = get_edge_sequence(node_set[root_child].parents[j], root_child);
				node_in = node_in + compute_coverage(kmer_hash, edge_in);
			}

			for (int j = 0; j < node_set[root_child].children.size(); j++) {
				const string& edge_out = get_edge_sequence(root_child, node_set[root_child].children[j]);
				node_out = node_out + compute_coverage(kmer_hash, edge_out);
			}


			TreeNode leaf_node(root_child, node_out-node_in);	

			
			for (int j = 0; j <node_set[root_child].cov_reads.size(); j++) {
				int mate_id = node_set[root_child].cov_reads[j].first;
				int rid;
				if (mate_id >= max_read_id)
					rid = mate_id - max_read_id;
				else
					rid = mate_id + max_read_id;
				if ( data_tag[rid] >= 0 && data_tag[rid] != root_child && (!node_set[root_child].is_child(data_tag[rid])) ) {
					leaf_node.add_candi(data_tag[rid], 1);
				}
			}
	
			int leaf_id = isotree.add_node(leaf_node);		
			isotree.add_child(root_id, leaf_id);
			pair <node_id, float> leaf_info;
			leaf_info.first = root_child;	
			leaf_info.second = node_out-node_in;	
			leaf.push_back(leaf_info);	
			tree_leaf_id.push_back(leaf_id);
		} else {
			for (int j = 0; j <node_set[root_child].cov_reads.size(); j++) {
				int rid = node_set[root_child].cov_reads[j].first; 
				data_tag[rid] = -1;
			}
		}

	}

	for (int i = 1; i < (int)node_order.size()-1; i++) {

		node_id x = node_order[node_order.size()-1-i];

/*
		if (node_set[x].node_coverage <= g_min_ratio_non_error ){
			for (int j = 0; j <node_set[x].cov_reads.size(); j++)
				data_tag[node_set[x].cov_reads[j].first] = -1;
			continue;
		}
*/
		dealed.push_back(x);

		float total_out = 0.0;
		vector<pair<node_id, float> > out_edge;
		bool has_t = false;

		for (int j = 0; j < node_set[x].children.size(); j++) {

			const string& edge = get_edge_sequence(x, node_set[x].children[j]);
			float edge_cov = compute_coverage(kmer_hash, edge);
			
			if (node_set[x].children[j] == node_order[0]) 
				has_t = true;			

			//if (node_set[node_set[x].children[j]].get_node_coverage() > g_min_ratio_non_error) {	//		 					
				total_out = total_out + edge_cov;								
				pair<node_id, float> edge_out_t;								
				edge_out_t.first = node_set[x].children[j];								
				edge_out_t.second = edge_cov;								
				out_edge.push_back(edge_out_t);				
			//} //


		}
		//if (total_out == 0.0)
			//has_t = true;

		float total_in = 0.0;
		vector<pair<node_id, float> > in_edge;
		vector<vector<int> > in_candi;
		vector<vector<float> > in_priority;
		for (int j = 0; j < leaf.size(); j++) {

			if (leaf[j].first == x) {  

				pair<node_id, float> in_edge_t;
				in_edge_t.first = j;
				in_edge_t.second = leaf[j].second;
				in_edge.push_back(in_edge_t);
				in_candi.push_back(isotree.treenode_set[tree_leaf_id[j]].candi);
				in_priority.push_back(isotree.treenode_set[tree_leaf_id[j]].priority);
				total_in = total_in + leaf[j].second;

			}

		}

		
		if (has_t) {
			pair<node_id, float> edge_out_t;								
			edge_out_t.first = node_order[0];								
			edge_out_t.second = total_in - total_out;								
			out_edge.push_back(edge_out_t);
			total_out = total_in;
		}
		
		if (total_out == 0.0 && total_in == 0.0) {

			if (node_set[x].sequence.length() >= g_min_transcript_length) {
			
				pair<string,float> transcript;			
				transcript.first = node_set[x].sequence;			
				transcript.second = node_set[x].node_coverage;			
				transcripts.push_back(transcript);
				vector<int> trans_node;
				trans_node.push_back(x);
				trans_node_vec.push_back(trans_node);
			}

			continue;

		}

		if ( total_out*(1-g_tolerance_value) > total_in*(1+g_tolerance_value)) {

			node_set[x].add_parent(node_order[node_order.size()-1]);
			node_set[node_order[node_order.size()-1]].add_child(x);
			TreeNode leaf_node(x, total_out-total_in);
			for (int it = 0; it > node_set[x].cov_reads.size();  it++) {
				int rid ;
				int mate_id = node_set[x].cov_reads[it].first;
				if (mate_id >= max_read_id)
					rid = mate_id - max_read_id;
				else
					rid = mate_id + max_read_id;
				if (data_tag[rid] >= 0 && find(dealed.begin(), dealed.end(), data_tag[rid]) != dealed.end() && (!node_set[x].is_child(data_tag[rid])) ) {
					leaf_node.add_candi(data_tag[rid], 1);
				}
			}
			int leaf_id = isotree.add_node(leaf_node);		
			isotree.add_child(root_id, leaf_id);
			pair <node_id, float> leaf_info;
			leaf_info.first = x;	
			leaf_info.second = total_out-total_in;	
			leaf.push_back(leaf_info);	
			tree_leaf_id.push_back(leaf_id);
			pair<node_id, float> in_edge_t;
			in_edge_t.first = leaf.size()-1;
			in_edge_t.second = total_out-total_in;
			in_edge.push_back(in_edge_t);
			in_candi.push_back(leaf_node.candi);
			in_priority.push_back(leaf_node.priority);
			total_in = total_out;

		}

		if ( total_out*(1+g_tolerance_value) < total_in*(1-g_tolerance_value)) {

			if (!has_t){

				node_set[x].add_child(node_order[0]);
				node_set[node_order[0]].add_parent(x);

			}
			pair <node_id, float> edge_out_t;
			edge_out_t.first = node_order[0];
			edge_out_t.second = total_in - total_out;
			out_edge.push_back(edge_out_t);
			total_out = total_in;

		}
		int out_sum = out_edge.size();
		int in_sum = in_edge.size();
		vector<float> adj(out_sum*in_sum, 0);
		solve_milp(out_edge, in_edge, adj, in_candi, in_priority);

		for (int j = 0; j < out_sum; j++) {

			//float re_sum = 0;

			//for (size_t k = 0; k < in_sum; k++)
				//re_sum += adj[k*out_sum+j];

			for (int k =0; k < in_sum; k++) {

				//adj[k*out_sum+j] =out_edge[j].second * adj[k*out_sum+j] / re_sum;

				if (adj[k*out_sum+j] > 0) {
				
					TreeNode c_node(out_edge[j].first, adj[k*out_sum+j]);
					c_node.candi = isotree.treenode_set[tree_leaf_id[in_edge[k].first]].candi;
					c_node.priority = isotree.treenode_set[tree_leaf_id[in_edge[k].first]].priority;
					c_node.dele_candi(out_edge[j].first);
					for (int it = 0; it > node_set[out_edge[j].first].cov_reads.size(); it++) {
						int rid;
						int mate_id = node_set[out_edge[j].first].cov_reads[it].first;
						if (mate_id >= max_read_id)
							rid = mate_id - max_read_id;
						else
							rid = mate_id + max_read_id;
						if (data_tag[rid] >= 0 && data_tag[rid] != out_edge[j].first &&  find(dealed.begin(), dealed.end(), data_tag[rid]) != dealed.end() && (!node_set[out_edge[j].first].is_child(data_tag[rid]))) {
							c_node.add_candi(data_tag[rid], 1);
						}
					}
				    int c_id = isotree.add_node(c_node);
				    isotree.add_child(tree_leaf_id[in_edge[k].first],c_id);
				    pair<int,float> leaf_x;
				    leaf_x.first = out_edge[j].first;
				    leaf_x.second = adj[k*out_sum+j];
				    leaf.push_back(leaf_x);
				    tree_leaf_id.push_back(c_id);
				}

			}

		}

		
		for (int j = 0 ; j < in_sum; j++) {

			for (int k = 0; k <in_sum-j-1; k++) {

				if (in_edge[k].first < in_edge[k+1].first) {

					pair<int, float> temp;
					temp = in_edge[k];
					in_edge[k] = in_edge[k+1];
					in_edge[k+1] = temp;
				
				}

			}

		}

		for(int j = 0; j < in_sum; j++) {

			vector<pair<node_id,float> >::iterator it = leaf.begin();
			leaf.erase(it+in_edge[j].first);
			vector<int>::iterator its = tree_leaf_id.begin();
			tree_leaf_id.erase(its+in_edge[j].first);

		}

	}
	

    if (isotree.total_node() > 1) {
	
		vector<float> w_trans;	
		vector<vector<int> > trans_vec;	
		isotree.dfs_trans(w_trans, trans_vec);
                
		for (int i = 0; i < w_trans.size(); i++) {
		
			string tran_seq;
		
			for (int j = 0; j < trans_vec[i].size(); j++) {

				if (trans_vec[i][j] == node_order[0] || trans_vec[i][j] == node_order[node_order.size()-1])
					continue;


				if (j > 0) {				
					tran_seq = tran_seq + node_set[trans_vec[i][j]].sequence;
					//cout << trans_vec[i][j] << ", "; ////
				}
				else {				
					tran_seq = node_set[trans_vec[i][j]].sequence;
					//cout << trans_vec[i][j] << ", "; ////
				}

		
			}
			if (w_trans[i] < g_min_trans_cov || tran_seq.length() < g_min_transcript_length) {
				for (int j = 0; j < trans_vec[i].size(); j++) {
					for (int it = 0; it < node_set[trans_vec[i][j]].cov_reads.size(); it++)
						data_tag[node_set[trans_vec[i][j]].cov_reads[it].first] = -1;	
				}
				continue;
			}
			pair<string, float> trans_x;		
			trans_x.first = tran_seq;		
			trans_x.second = w_trans[i];	
			transcripts.push_back(trans_x);
			trans_node_vec.push_back(trans_vec[i]);
	
		}
		
	}

	//check_transcripts(kmer_hash, transcripts, trans_node_vec);
	check_transcripts(kmer_hash, transcripts);
}


void SplicingGraph::topological_sort() {
	
	Node s,t;
	s.sequence = "start";
	t.sequence = "end";
	int s_id = add_node(s);
	int t_id = add_node(t);
	vector<int> node_color;

	for (int i = 0; i < node_sum; i++) {

		if (node_set[i].parents.size() == 0 && i != s_id && i != t_id){//

			node_set[i].add_parent(s_id);
			node_set[s_id].add_child(i);

		} 

		if (node_set[i].children.size() == 0 && i != s_id && i != t_id) {//

			node_set[i].add_child(t_id);
			node_set[t_id].add_parent(i);

		}

		node_color.push_back(0);

	}

	dfs_visit(s_id, node_color);

}


void SplicingGraph::dfs_visit(node_id i, vector<int>& node_color) {

	node_color[i] = 1;

	if (node_set[i].children.size() == 0) 
		node_order.push_back(i);
	else {

		for (int j = 0 ; j < node_set[i].children.size(); j++) {

			if (node_color[node_set[i].children[j]] == 0)
				dfs_visit(node_set[i].children[j], node_color);

		}

		node_order.push_back(i);

	}

}



void SplicingGraph::describe_graph(KmerHash& kmer_hash, size_t splicing_graph_id) {

	vector<float> edge_m(node_sum*node_sum,0);
	stringstream top_node;
	top_node <<"topogical_node" << splicing_graph_id;
	string top_name ="G:\\Jane\\" + top_node.str() + ".txt";
	fstream topgical_file;
	topgical_file.open(top_name.c_str(), fstream::out);

	if (!topgical_file.is_open()) {

		cout << " File " << top_name.c_str() << " can't be opened!" << endl;
		exit(1);

	}

	for (int i = 0; i < node_order.size(); i++) 
		topgical_file << node_order[node_order.size()-1-i] << endl;
	
	topgical_file.close();

	stringstream node_seq;
	node_seq <<"node_sequence" << splicing_graph_id;
	string node_name = "G:\\Jane\\" + node_seq.str() + ".txt";
	fstream node_file;
	node_file.open(node_name.c_str(), fstream::out);

	if (!node_file.is_open()) {

		cout << " File " << node_name.c_str() << " can't be opened!" << endl;
		exit(1);

	}

	for (int i = 0; i < node_order.size(); i++) {

		node_id x = node_order[node_order.size()-1-i];

		node_file << "node id = " << x << "\tlength = " << node_set[x].sequence.length();

		string sequence = node_set[x].sequence;
		float node_cov;

		if ((int)sequence.length() < g_read_length) //the node sequence is too short
			node_cov = 0.0;
		else
			node_cov = compute_coverage(kmer_hash, sequence);

		edge_m[x*node_sum + x]= node_cov;

		//if (sequence.length() > g_read_length && node_cov == 0.0)
			//recover_kmers(sequence);

		node_file << "\tcov = " << setiosflags(ios::fixed) << setprecision(2) << node_cov 
			<< "\tsequence : " << endl << node_set[x].sequence << endl;

	}
	node_file.close();

	stringstream edge_matrix;
	edge_matrix <<"edge_matrix" << splicing_graph_id;
	string edge_name = "G:\\Jane\\" + edge_matrix.str() + ".txt";
	fstream edge_file;
	edge_file.open(edge_name.c_str(), fstream::out);

	if (!edge_file.is_open()) {

		cout << " File " << edge_name.c_str() << " can't be opened!" << endl;
		exit(1);

	}

	for (int i = 0; i < node_sum; i++) {

		for (int j = 0; j < node_set[i].children.size(); j++) {
			
			string edge = get_edge_sequence(i, node_set[i].children[j]);
			edge_m[i*node_sum+node_set[i].children[j]] = compute_coverage(kmer_hash, edge);

		}

	}

	for (int i = 0; i < edge_m.size(); i++) {

		edge_file << edge_m[i] << " ";

		if ((i+1)%node_sum == 0 && i > 0)
			edge_file << endl;

	}

	edge_file.close();

}


int SplicingGraph::get_total_amount_of_kmers() {

	int kmer_count = 0;

	for (unsigned int i = 0; i < node_sum; i++)
		kmer_count += node_set[i].sequence.length();

	return (kmer_count - g_kmer_length);

}





void SplicingGraph::solve_milp(vector<pair<node_id, float> >& out_edge, vector<pair<node_id, float> >& in_edge, vector<float>& adj, vector<vector<int> >& in_candi, vector<vector <float> >& in_priority) {
		
	int out_sum = out_edge.size();
	int in_sum = in_edge.size();
	if (out_sum == 1 || in_sum == 1) {
		if (out_sum == 1 && in_sum == 1)			
			adj[0] = out_edge[0].second;
		else {
			if (out_sum == 1) {
				for (int j = 0; j < in_sum; j++)
					adj[j] = in_edge[j].second;
			}
			if (in_sum == 1) {
				for (int j = 0; j < out_sum; j++)
					adj[j] = out_edge[j].second;
			}
		}

	} else {
		bool assinged = false;
		int in_t = 0;
		int out_t = 0;
		vector<float> in_relax(in_sum);
		vector<float> out_relax(out_sum);
		vector<int> in_tag(in_sum);
		vector<int> out_tag(out_sum);
		//新方法：
		for (int i = 0; i < out_sum; i++) {
			float max_candi = 0;
			int candi_id = -1;
			for (int j = 0; j < in_sum; j++) {

				for (int k = 0; k < in_candi[j].size(); k++) {
					if (in_candi[j][k] == out_edge[i].first) {
						if (in_priority[j][k] > max_candi && in_tag[j] == 0) {
							max_candi = in_priority[j][k];
							candi_id = j;
						}
						break;
					}
				}
			}

			if (max_candi > 0) {

				if (out_edge[i].second > in_edge[candi_id].second)					
					adj[out_sum*candi_id+i] = in_edge[candi_id].second;
				else
					adj[out_sum*candi_id+i] = out_edge[i].second;

				in_tag[candi_id] =1;
				out_tag[i] = 1;
			}

		}
		
		for (int i = 0; i < out_sum; i++) {

			if (out_tag[i] == 0) {
				bool assign_tag = false;
				for (int j = 0; j < in_sum; j++) {
					if (in_tag[j] == 0 && assign_tag == false) {
						if (in_edge[j].second*(1-g_tolerance_value) <= out_edge[i].second*(1+g_tolerance_value) && in_edge[j].second*(1+g_tolerance_value)>=out_edge[i].second*(1-g_tolerance_value)) {
							adj[out_sum*j+i] = in_edge[j].second;
							in_tag[j] = 1;
							out_tag[i] = 1;
							in_relax[j] = in_edge[j].second*g_tolerance_value;
							out_relax[i] = out_edge[i].second*(1+g_tolerance_value) - in_edge[j].second;
							assign_tag = true;
							break;
						} else {
							if (in_edge[j].second > out_edge[i].second) {
								adj[out_sum*j+i] = out_edge[i].second;
								out_tag[i] = 1;
								assign_tag = true;
								out_relax[i] = out_edge[i].second*g_tolerance_value;
				                in_edge[j].second = in_edge[j].second - out_edge[i].second;
								break;
							} else {
								adj[out_sum*j+i] = in_edge[j].second;
								in_tag[j] = 1;
								in_relax[j] = in_edge[j].second*g_tolerance_value;
								out_edge[i].second = out_edge[i].second - in_edge[j].second;
							}
						}
					}
				}
				if (assign_tag == false) {
					for (int j = 0; j < in_sum; j++) {
						if (in_relax[j] <= out_edge[i].second*(1+g_tolerance_value) && in_relax[j]>=out_edge[i].second*(1-g_tolerance_value)) {
							adj[out_sum*j+i] = in_relax[j];							
							out_tag[i] = 1;
							out_relax[i] = out_edge[i].second*(1+g_tolerance_value) - in_relax[j];
							in_relax[j] = 0;
							assign_tag = true;
							break;
						} else {
							if (in_relax[j] > out_edge[i].second) {
								adj[out_sum*j+i] = out_edge[i].second;
								out_tag[i] = 1;
								assign_tag = true;
								out_relax[i] = out_edge[i].second*g_tolerance_value;
				                in_relax[j] = in_relax[j] - out_edge[i].second;
								break;
							} else {
								adj[out_sum*j+i] = in_relax[j];
								out_edge[i].second = out_edge[i].second - in_relax[j];
								in_relax[j] = 0;
							}
						}				
					}
				}

			}

		}


		
		for (int i = 0; i < in_sum; i++) {

			if (in_tag[i] == 0) {
				bool assign_tag = false;
				for (int j = 0; j < out_sum; j++) {
					if (out_tag[j] == 0 && assign_tag == false) {
						if (out_edge[j].second*(1-g_tolerance_value) <= in_edge[i].second*(1+g_tolerance_value) && out_edge[j].second*(1+g_tolerance_value)>=in_edge[i].second*(1-g_tolerance_value)) {
							adj[out_sum*i+j] = out_edge[j].second;
							out_tag[j] = 1;
							in_tag[i] = 1;
							out_relax[j] = out_edge[j].second*g_tolerance_value;
							in_relax[i] = in_edge[i].second*(1+g_tolerance_value) - out_edge[j].second;
							assign_tag = true;
							break;
						} else {
							if (out_edge[j].second > in_edge[i].second) {
								adj[out_sum*i+j] = in_edge[i].second;
								in_tag[i] = 1;
								assign_tag = true;
								in_relax[i] = in_edge[i].second*g_tolerance_value;
				                out_edge[j].second = out_edge[j].second - in_edge[i].second;
								break;
							} else {
								adj[out_sum*i+j] = out_edge[j].second;
								out_tag[j] = 1;
								out_relax[j] = out_edge[j].second*g_tolerance_value;
								in_edge[i].second = in_edge[i].second - out_edge[j].second;
							}
						}
					}
				}
				if (assign_tag == false) {
					for (int j = 0; j < out_sum; j++) {
						if (out_relax[j] <= in_edge[i].second*(1+g_tolerance_value) && out_relax[j]>=in_edge[i].second*(1-g_tolerance_value)) {
							adj[out_sum*i+j] = out_relax[j];
							in_tag[i] = 1;
							in_relax[i] = in_edge[i].second*(1+g_tolerance_value) - out_relax[j];
							out_relax[j] = 0;
							assign_tag = true;
							break;
						} else {
							if (out_relax[j] > in_edge[i].second) {
								adj[out_sum*i+j] = in_edge[i].second;
								in_tag[i] = 1;
								assign_tag = true;
								in_relax[i] = in_edge[i].second*g_tolerance_value;
				                out_relax[j] = out_relax[j] - in_edge[i].second;
								break;
							} else {
								adj[out_sum*i+j] = out_relax[j];
								in_edge[i].second = in_edge[i].second - out_relax[j];
								out_relax[j] = 0;
							}
						}				
					}
				}

			}

		}

		vector<float> sum_vec;
		for (int i=0; i< out_sum; i++)
			sum_vec.push_back(out_edge[i].second);

		for (int i = 0; i < out_sum; i++) {
			float re_sum = 0.0;
			for (int j = 0; j < in_sum; j++) 
				re_sum = re_sum + adj[j*out_sum+i];
			for (int j = 0; j < in_sum; j++) {
				adj[j*out_sum+i] = sum_vec[i]*adj[j*out_sum+i] / re_sum;
			}
		}
	}
		
}


void SplicingGraph::recover_reads() {

	for (int i = 0; i < node_sum; i++) {
		
		for (int it = 0; it < node_set[i].cov_reads.size(); it++) {
			data_tag[node_set[i].cov_reads[it].first] = -1;
		}
	}



}

void SplicingGraph::check_transcripts(KmerHash& kmer_hash, vector<pair<string,float> >& transcripts, vector<vector<int> >& trans_vec) {
	cout << "check_transcripts.." << endl;
	vector<int> add_trans;
	vector<int> set_reads_vec;
	vector<int> recover_reads_vec;
	set<int> used_reads;
	for (int i = 0; i < transcripts.size(); i++) {

		string transcript = transcripts[i].first;
		bool is_add = true;
		if (transcripts[i].second < g_min_trans_cov || transcript.length() < g_min_transcript_length) {
			for (int j = 0; j < trans_vec[i].size(); j++) {
				for (int it = 0; it < node_set[trans_vec[i][j]].cov_reads.size(); it++) {
					int rid = node_set[trans_vec[i][j]].cov_reads[it].first;
					if (data_tag[rid] != -2)
						data_tag[rid] = -1;
				}
			}
			continue;
		}
		
		if (g_is_paired_end) {
			list<int> reads_in_trans1;
			vector<int> reads_in_trans;
			int pair_end_sum = 0;			
			for (int j = 0; j < trans_vec[i].size(); j++) {
				for (int it = 0; it < node_set[trans_vec[i][j]].cov_reads.size(); it++) {
					int rid = node_set[trans_vec[i][j]].cov_reads[it].first;
					if (data_tag[rid] == -2)//
						continue;//
					data_tag[rid] = -4;
					reads_in_trans.push_back(rid);
				}
			}
			vector<int> edge_reads;
			for (int j = 0; j < trans_vec[i].size()-1; j++) {				
				get_edge_reads(kmer_hash, trans_vec[i][j], trans_vec[i][j+1], edge_reads);
			}
			for (int j = 0; j < edge_reads.size(); j++) {
				if (data_tag[edge_reads[j]] == -2)//
					continue;//
				data_tag[edge_reads[j]] = -4;
				reads_in_trans.push_back(edge_reads[j]);
			}
			vector<int> maybe_used;
			for (int j = 0; j < reads_in_trans.size(); j++) {

				int read_id = reads_in_trans[j];
				int mate_id;
				if (read_id >= max_read_id)
					mate_id = read_id - max_read_id;
				else
					mate_id = read_id + max_read_id;

				if (data_tag[read_id] == -4 && data_tag[mate_id] == -4) {
					maybe_used.push_back(read_id);
					maybe_used.push_back(mate_id);
					data_tag[mate_id] = -1;
					pair_end_sum = pair_end_sum + 2;					
				}

				data_tag[read_id] = -1;
			}

			transcripts[i].second = static_cast<float>(pair_end_sum)/ (transcripts[i].first.length()-g_read_length+1);
			float pair_ratio = static_cast<float>(pair_end_sum)/reads_in_trans.size();

			if (pair_ratio < g_min_pair_ratio) {
				is_add = false;
			} else {
				for (int j = 0; j < maybe_used.size(); j++)
					used_reads.insert(maybe_used[j]);
			}
		}

		if (is_add)
			add_trans.push_back(i);

	}
	

	vector<pair<string,float> > temp_trans;
	for (int i = 0; i < add_trans.size(); i++) {

		int tid = add_trans[i];
		temp_trans.push_back(transcripts[tid]);
	}

	set<int>::iterator it;
	for (it = used_reads.begin(); it != used_reads.end(); it++)
		data_tag[*it] = -2;

	transcripts.clear();
	for (int i = 0; i < temp_trans.size(); i++)
		transcripts.push_back(temp_trans[i]);

}


int SplicingGraph::get_read_id(node_id p, int pos){


	for (int i = pos; i < node_set[p].sequence.length()-g_read_length; i++) {

		string read = node_set[p].sequence.substr(i, g_read_length);
		for(int it = 0; it < node_set[p].cov_reads.size(); it++){
			int rid = node_set[p].cov_reads[it].first;
			if (data[rid] == read)
				return rid;
		}
	}

	return -1;
}


void SplicingGraph::get_edge_reads(KmerHash & kmer_hash, node_id source, node_id target, vector<int>& edge_reads) {

	if (source < 0 || source >= node_sum)
		return;
	if (target < 0 || target >= node_sum)
		return;
	string sequence = get_edge_sequence(source, target);
	for (int i = 0; i <= sequence.length()-g_read_length; i++) {
		string real_read = sequence.substr(i,g_read_length);
		kmer_int_type a_kmer = kmer_to_int(sequence.substr(i, g_kmer_length)); 
		kmer_int_type b_kmer = kmer_to_int(sequence.substr(i+g_read_length-g_kmer_length, g_kmer_length));
		vector<int> a_readset = kmer_hash[a_kmer];
		vector<int> b_readset = kmer_hash[b_kmer];
		vector<int>::iterator a = a_readset.begin();
		vector<int>::iterator b = b_readset.begin();
		while (a != a_readset.end() && b != b_readset.end()) {

			if (*a == *b) {
				if (data_tag[*a] != -2) {
					string read = data[*a];				
					if (read == real_read)				
						edge_reads.push_back(*a);
					else {
						string read_rev = revcomp(read);
				    	if (g_is_paired_end && g_double_stranded_mode && read_rev == real_read)
					    	edge_reads.push_back(*a);
					}
				}
				a++;
				b++;

			} else {

				if (*a > *b)
					b++;
				else
					a++;

			}

		}

	}


}



void SplicingGraph::set_reads_tag(KmerHash& kmer_hash, string sequence, int tag) {

	if (sequence.length() < g_read_length)
		return;

	for (int i = 0; i <= sequence.length() - g_read_length; i++) {
		string real_read = sequence.substr(i,g_read_length);
		kmer_int_type a_kmer = kmer_to_int(sequence.substr(i, g_kmer_length)); 
		kmer_int_type b_kmer = kmer_to_int(sequence.substr(i+g_read_length-g_kmer_length, g_kmer_length));
		vector<int> a_readset = kmer_hash[a_kmer];
		vector<int> b_readset = kmer_hash[b_kmer];
		vector<int>::iterator a = a_readset.begin();
		vector<int>::iterator b = b_readset.begin();
		while (a != a_readset.end() && b != b_readset.end()) {

			if (*a == *b) {
				if (data_tag[*a] != -2) {
					string read = data[*a];				
					if (read == real_read)				
						data_tag[*a] = tag;
					else {
						string read_rev = revcomp(read);
						if (g_is_paired_end && g_double_stranded_mode && read_rev == real_read)
							data_tag[*a] = tag;
					}
				}
				a++;
				b++;

			} else {

				if (*a > *b)
					b++;
				else
					a++;
	
			}

		}
	}
}


void SplicingGraph::set_reads_tag(KmerHash& kmer_hash, string sequence, list<int>& map_reads, int tag) {

	if (sequence.length() < g_read_length)
		return;

	for (int i = 0; i <= sequence.length() - g_read_length; i++) {
		string real_read = sequence.substr(i,g_read_length);
		kmer_int_type a_kmer = kmer_to_int(sequence.substr(i, g_kmer_length)); 
		kmer_int_type b_kmer = kmer_to_int(sequence.substr(i+g_read_length-g_kmer_length, g_kmer_length));
		vector<int> a_readset = kmer_hash[a_kmer];
		vector<int> b_readset = kmer_hash[b_kmer];
		vector<int>::iterator a = a_readset.begin();
		vector<int>::iterator b = b_readset.begin();
		while (a != a_readset.end() && b != b_readset.end()) {

			if (*a == *b) {
				if (data_tag[*a] != -2) {
					string read = data[*a];				
					if (read == real_read) {				
						data_tag[*a] = tag;
						map_reads.push_back(*a);
					}
					else {
						string read_rev = revcomp(read);
						if (g_is_paired_end && g_double_stranded_mode && read_rev == real_read){
							data_tag[*a] = tag;
							map_reads.push_back(*a);
						}
					}
				}
				a++;
				b++;

			} else {

				if (*a > *b)
					b++;
				else
					a++;
	
			}

		}
	}
}




void SplicingGraph::map_reads_to_sequence(KmerHash& kmer_hash, string sequence, vector<pair<int, int> >& reads_pos) {

	if (sequence.length() < g_read_length)
		return;
	
	for (int i = 0; i <= sequence.length() - g_read_length; i++) {
		string real_read = sequence.substr(i,g_read_length);
		kmer_int_type a_kmer = kmer_to_int(sequence.substr(i, g_kmer_length)); 
		kmer_int_type b_kmer = kmer_to_int(sequence.substr(i+g_read_length-g_kmer_length, g_kmer_length));
		vector<int> a_readset = kmer_hash[a_kmer];
		vector<int> b_readset = kmer_hash[b_kmer];
		vector<int>::iterator a = a_readset.begin();
		vector<int>::iterator b = b_readset.begin();
		while (a != a_readset.end() && b != b_readset.end()) {

			if (*a == *b) {
				if (data_tag[*a] != -2) {
					string read = data[*a];				
					if (read == real_read){				
						pair<int, int> r_p;
						r_p.first = *a;
						r_p.second = i;
						data_tag[*a] = i;
						reads_pos.push_back(r_p);
					} else {
						string read_rev = revcomp(read);
						if (g_is_paired_end && g_double_stranded_mode && read_rev == real_read) {
							pair<int, int> r_p;
							r_p.first = *a;
							r_p.second = i;
							data_tag[*a] = i;
							reads_pos.push_back(r_p);
						}
					}
				}
				a++;
				b++;

			} else {

				if (*a > *b)
					b++;
				else
					a++;
	
			}

		}
	}

}


void SplicingGraph::check_transcripts(KmerHash& kmer_hash, vector<pair<string,float> >& transcripts) {

	cout << "check_transcripts.." << endl;

	vector<int> add_trans;
	vector<int> set_reads_vec;
	vector<int> recover_reads_vec;
	set<int> used_reads;
///*
	string max_transcript = "N";
	int max_tran_id = -1;
	for (int i = 0; i < transcripts.size(); i++) 
		if (transcripts[i].first.length() > max_transcript.length()) {
			max_transcript = transcripts[i].first;
			max_tran_id = i;
		}

//*/
	for (int i = 0; i < transcripts.size(); i++) {
		string transcript = transcripts[i].first;
		vector<pair<int, int> > reads_pos;
		map_reads_to_sequence(kmer_hash, transcript, reads_pos);
		bool is_add = true;
///*
		if (is_seq_similar(max_transcript, transcript, 'F',0.05) && i != max_tran_id)
			is_add = false;
		if (is_seq_similar(max_transcript, transcript,'R', 0.05) && i != max_tran_id)
			is_add = false;
//*/
cout << transcripts[i].second << endl;
		if ((!is_add) || transcripts[i].second < g_min_trans_cov || transcript.length() < g_min_transcript_length) {

			for (int j = 0; j < reads_pos.size(); j++)
				data_tag[reads_pos[j].first] = -1;
			continue;
		}
		
		if (g_is_paired_end) {
			
			float pair_end_sum = 0.0;
/*			
			for (int j = 0; j < reads_pos.size(); j++) {
				if (data_tag[reads_pos[j].first] == -2)
					continue;
				data_tag[reads_pos[j].first] = -4;
			}
*/			
			vector<pair<int, int> > maybe_used;
			for (int j = 0; j < reads_pos.size(); j++) {

				int read_id = reads_pos[j].first;
				int mate_id;
				if (read_id >= max_read_id)
					mate_id = read_id - max_read_id;
				else
					mate_id = read_id + max_read_id;

				if (data_tag[read_id] >= 0 && data_tag[mate_id] >= 0) {
					maybe_used.push_back(reads_pos[j]);
					pair_end_sum = pair_end_sum + 1.0;
				} else {
					data_tag[read_id] = -1;
				}
/*
				if (read_id < mate_id && data_tag[read_id] <= data_tag[mate_id]) {
					maybe_used.push_back(reads_pos[j]);
					pair_end_sum = pair_end_sum + 1.0;
				} else {
					if (read_id > mate_id && data_tag[read_id] >= data_tag[mate_id] && data_tag[mate_id] >= 0) {
						maybe_used.push_back(reads_pos[j]);
						pair_end_sum = pair_end_sum + 1.0;
					}else {
						data_tag[read_id] = -1;
					}
				}
*/
/*
				if (data_tag[read_id] == -4 && data_tag[mate_id] == -4) {

					maybe_used.push_back(reads_pos[j]);
					pair_end_sum = pair_end_sum + 1.0;					
				} else 
					data_tag[read_id] = -1;
*/
			}
cout << pair_end_sum << " " << reads_pos.size() << endl;

			if (pair_end_sum > 10 && pair_end_sum / reads_pos.size() > 0.4) {
				int total_cov_len = 0;
				int current_cov_len = g_read_length;
				int current_cov_start = maybe_used[0].second;
				bool add_tag = true;
				for (int j = 0; j < maybe_used.size()-1; j++) {
					if (maybe_used[j].second + g_read_length >= maybe_used[j+1].second) {
						add_tag = true;
						current_cov_len = maybe_used[j+1].second + g_read_length - current_cov_start;
					} else {
						add_tag = false;
						total_cov_len = total_cov_len + current_cov_len;
						current_cov_len = g_read_length;
						current_cov_start = maybe_used[j+1].second;
					}
				}
				if (add_tag = true)
					total_cov_len = total_cov_len + current_cov_len;
				if (total_cov_len > 0.75*transcript.length()){
					vector<int> trans_cov(transcript.length(), 0);
					for (int j = 0; j < maybe_used.size(); j++) {
						for (int k = maybe_used[j].second; k < maybe_used[j].second + g_read_length; k++)
							trans_cov[k] = trans_cov[k]+1;
					}
					vector<int> xxx = trans_cov;
					sort(trans_cov.begin(), trans_cov.end());
					int check_id = trans_cov.size() * 0.4;
					int min_cov = xxx[check_id];
					if (min_cov < transcripts[i].second)
						min_cov = transcripts[i].second;				
					//int min_cov = transcripts[i].second;
					vector<int> trans_cov1 = trans_cov;
					int current_check_id = maybe_used[0].first;
					for (int j = 1; j < maybe_used.size(); j++){
						if (data[current_check_id] == data[maybe_used[j].first])
							continue;
						else
							current_check_id = maybe_used[j].first;
						int support_times = 0;
						for (int k = maybe_used[j].second; k < maybe_used[j].second + g_read_length; k++) {
							if (trans_cov[k] - trans_cov1[k] <= min_cov || trans_cov[k] <= min_cov) {
								support_times = support_times + 1;
							}
						}
						if (support_times >= g_read_length-1) {
							used_reads.insert(maybe_used[j].first);
							for (int t = maybe_used[j].second; t < maybe_used[j].second + g_read_length; t++)
								trans_cov1[t] = trans_cov1[t] - 1;
							int mate_id;
							if (maybe_used[j].first >= max_read_id)
								mate_id = maybe_used[j].first-max_read_id;
							else
								mate_id =maybe_used[j].first+max_read_id;
							used_reads.insert(mate_id);
							for (int t = data_tag[mate_id]; t < data_tag[mate_id] + g_read_length; t++){
								if (t < 0)
									break;
								trans_cov1[t] = trans_cov1[t] - 1;
							}
						}
					}

					for (int j = 0; j < maybe_used.size(); j++)
						data_tag[maybe_used[j].first] = -1;
/*	
					for (int j = 0; j < maybe_used.size(); j++) 
						used_reads.insert(maybe_used[j].first);
*/
					transcripts[i].second = pair_end_sum * g_read_length / transcript.length();
				} else {
					is_add = false;
					for (int j = 0; j < maybe_used.size(); j++)
						data_tag[maybe_used[j].first] = -1;
				}
			} else { 
				is_add = false;
				for (int j = 0; j < maybe_used.size(); j++)
					data_tag[maybe_used[j].first] = -1;
			}
					
		}
		if (is_add) {
			add_trans.push_back(i);
/*
			set<int>::iterator it;
			for (it = used_reads.begin(); it != used_reads.end(); it++)
			data_tag[*it] = -2;
			used_reads.clear();
*/
		}

	}
	

	vector<pair<string,float> > temp_trans;
	pair<string,float> max_trans;
	max_trans.first = " ";
	max_trans.second = 0.0;
	for (int i = 0; i < add_trans.size(); i++) {

		int tid = add_trans[i];
		if (transcripts[tid].first.length() > max_trans.first.length())
			max_trans = transcripts[tid];
		temp_trans.push_back(transcripts[tid]);
	}
//*
	set<int>::iterator it;
	for (it = used_reads.begin(); it != used_reads.end(); it++)
		data_tag[*it] = -2;
//*/
	transcripts.clear();	
	transcripts = temp_trans;

}


void SplicingGraph::refine_transcript(vector<pair<int, int> > maybe_used, set<int>& used_reads, int& start_pos, int& trans_len) {

	start_pos = 0;
	trans_len = 0;
	if (maybe_used.size() < 10){
        for (int j = 0; j < maybe_used.size(); j++) 
            data_tag[maybe_used[j].first] = -1;
		return;
    }
	int max_cov_len = g_read_length;
	int max_cov_start = maybe_used[0].second;
	int current_cov_start = maybe_used[0].second;
	int current_cov_len = g_read_length;
	int refine_length = maybe_used[maybe_used.size()-1].second - maybe_used[0].second + g_read_length;

    if (refine_length < g_min_transcript_length) {
        for (int j = 0; j < maybe_used.size(); j++) 
            data_tag[maybe_used[j].first] = -1;
        return;
    }

	for (int j = 0; j < maybe_used.size()-1; j++) {
					//cout << maybe_used[j].second <<",";
		if (maybe_used[j].second + g_read_length >= maybe_used[j+1].second) {
			current_cov_len = maybe_used[j+1].second + g_read_length - current_cov_start;
		} else {
			if (current_cov_len > max_cov_len){
				max_cov_len = current_cov_len;
				max_cov_start = current_cov_start;
			}
			current_cov_len = g_read_length;
			current_cov_start = maybe_used[j+1].second;
		}
	}
	if (current_cov_len > max_cov_len){
		max_cov_len = current_cov_len;
		max_cov_start = current_cov_start;
	}

	if (max_cov_len == refine_length) {
		//if (max_cov_len == transcript.length()){
		for (int k = 0; k < maybe_used.size(); k++) 
			used_reads.insert(maybe_used[k].first);
		//} else {
			//for (int k = 0; k < maybe_used.size(); k++) 
				//data_tag[maybe_used[k].first] = -1;
		//}
		start_pos = max_cov_start;
		trans_len = max_cov_len;
		
	} else {
		vector<pair<int, int> > new_used;
		for (int k = 0; k < maybe_used.size(); k++) {
			if (maybe_used[k].second < max_cov_start || maybe_used[k].second > max_cov_start+max_cov_len-g_read_length) {				
				data_tag[maybe_used[k].first] = -1;
			} else {
				if (data_tag[maybe_used[k].first] == -4)
					new_used.push_back(maybe_used[k]);
			}
		}
		vector<pair<int, int> > new_used1;
		for (int k = 0; k < new_used.size(); k++) {
			if (data_tag[new_used[k].first] == -4) {
				int mate_id;
				if (new_used[k].first >= max_read_id)
					mate_id = new_used[k].first - max_read_id;
				else
					mate_id = new_used[k].first + max_read_id;
				if (data_tag[mate_id] == -4)
					new_used1.push_back(new_used[k]);
                		else
                    			data_tag[new_used[k].first] = -1;
			}
		}
		refine_transcript(new_used1, used_reads, start_pos, trans_len);
	}

}


void SplicingGraph::set_reads_pos_in_node(KmerHash& kmer_hash, node_id p) {

	if (p < 0 || p >= node_sum )
		return;
	if (node_set[p].sequence.length() < g_read_length)
		return;

	string sequence = node_set[p].sequence;

	if (node_set[p].cov_reads.size() > 0) {
		for (int i = 0; i < node_set[p].cov_reads.size(); i++)
			data_tag[node_set[p].cov_reads[i].first] = -1;
		node_set[p].cov_reads.clear();
	}

	for (int i = 0; i <= sequence.length() - g_read_length; i++) {
		string real_read = sequence.substr(i,g_read_length);
		kmer_int_type a_kmer = kmer_to_int(sequence.substr(i, g_kmer_length)); 
		kmer_int_type b_kmer = kmer_to_int(sequence.substr(i+g_read_length-g_kmer_length, g_kmer_length));
		vector<int> a_readset = kmer_hash[a_kmer];
		vector<int> b_readset = kmer_hash[b_kmer];
		vector<int>::iterator a = a_readset.begin();
		vector<int>::iterator b = b_readset.begin();
		while (a != a_readset.end() && b != b_readset.end()) {

			if (*a == *b) {
				if (data_tag[*a] != -2) {
					string read = data[*a];				
					if (read == real_read){				
						pair<int, int> r_p;
						r_p.first = *a;
						r_p.second = i;
						node_set[p].cov_reads.push_back(r_p);
						data_tag[*a] =p;
					} else {
						string read_rev = revcomp(read);
						if (g_is_paired_end && g_double_stranded_mode && read_rev == real_read) {
							pair<int, int> r_p;
							r_p.first = *a;
							r_p.second = i;
							node_set[p].cov_reads.push_back(r_p);
							data_tag[*a] = p;
						}
					}
				}
				a++;
				b++;

			} else {

				if (*a > *b)
					b++;
				else
					a++;
	
			}

		}
	}


}



void SplicingGraph::branch_extend_by_coverage(KmerHash& kmer_hash) {

	forward_extend_by_coverage(kmer_hash, 0);

	int current_node_sum = node_sum;

	for (int i = 0; i < current_node_sum; i++)
		reverse_extend_by_coverage(kmer_hash, i);

}


void SplicingGraph::forward_extend_by_coverage(KmerHash& kmer_hash, node_id p) {

	if (p < 0 || p >= node_sum)
		return;
cout << "forward check " << p << endl;
	while (node_set[p].forward_check_pos.size() > 0) {

		bool break_while = false;
		pair<int, int> check_pos = node_set[p].forward_check_pos.front();
		node_set[p].forward_check_pos.pop_front();
cout << node_set[p].sequence.length() << " " << check_pos.first << " " << check_pos.second << endl;
//cout << check_pos.first << endl;
//cout << check_pos.second << endl;
//cout << node_set[p].sequence << endl;
		for (int i = check_pos.first; i < check_pos.second; i++) {
//cout << i << endl;
			string seed = node_set[p].sequence.substr(i, g_read_length);
			list<int> map_reads;
			pair<int, int> right_id;
			right_id.first = -1;
			right_id.second = -1;
			string forward_seq = forward_extend(kmer_hash, seed, map_reads, right_id);
			list<int>::iterator it;//
			for(it = map_reads.begin(); it != map_reads.end(); it++)//
				data_tag[*it] = -1;//

			if (right_id.first >= 0 && data_tag[right_id.first] >= node_sum) {
				right_id.first = -1;
				right_id.second = -1;
			}
			string::size_type start2;
			int right_node = -1;
			if (right_id.first >= 0) {
				right_node = data_tag[right_id.first];
				if (!pair_support(map_reads, right_node)) {
					right_node = -1;
					right_id.first = -1;
					right_id.second = -1;
				} 
			}
			if ((!pair_support(map_reads,p)) || forward_seq.length() <= g_read_length + 10) {
				//list<int>::iterator it;
				//for(it = map_reads.begin(); it != map_reads.end(); it++)
					//data_tag[*it] = -1;
				continue;
			}

			if (right_id.first == -1 && forward_seq.length() < g_min_exon_length + g_read_length) {
				//list<int>::iterator it;
				//for(it = map_reads.begin(); it != map_reads.end(); it++)
					//data_tag[*it] = -1;
				continue;
			}

			if (right_id.first == -1) {
//cout <<"-1" << endl;					
				//if (!is_seq_similar(forward_seq.substr(g_read_length), node_set[p].sequence.substr(i + g_read_length))) {

//cout << "similar" << endl;
//cout << node_set[p].sequence.length() << endl;
//cout << forward_seq.length() << endl;

					Node node2 = node_set[p].sequence.substr(0, i + g_read_length);
					int n2 = add_node(node2);
					node_set[p].sequence = node_set[p].sequence.substr(i+g_read_length);
					for (int j = 0; j < node_set[p].parents.size(); j++) {
						int p_i = node_set[p].parents[j];
						node_set[p_i].delete_child(p);
						node_set[p_i].add_child(n2);
					}
					node_set[n2].parents = node_set[p].parents;
					node_set[p].parents.clear();
					node_set[p].add_parent(n2);
					node_set[n2].add_child(p);
					if (node_set[p].sequence.length() < g_min_exon_length && node_set[p].children.size() == 0) {

//cout << "if " << endl;
//cout << i << endl;
						//set_reads_tag(kmer_hash, node_set[p].sequence, -1);
						node_set[p].sequence = forward_seq.substr(g_read_length);
						set_reads_pos_in_node(kmer_hash, p);
						set_check_pos(p);
						forward_extend_by_coverage(kmer_hash, p);
						break_while = true;
						break;
					} else {

//cout << "else" << endl;
//cout << i << endl;
						Node node1;
						node1.sequence = forward_seq.substr(g_read_length);
						int n1 = add_node(node1);
						node_set[n2].add_child(n1);
						node_set[n1].add_parent(n2);
						set_reads_pos_in_node(kmer_hash, p);
						set_reads_pos_in_node(kmer_hash, n1);
						set_reads_pos_in_node(kmer_hash, n2);
						set_check_pos(p);
						set_check_pos(n1);
						set_check_pos(n2);
						forward_extend_by_coverage(kmer_hash, p);
						//forward_extend_by_coverage(kmer_hash, n1);//
						break_while = true;
						break;
					}
						
				/*} else {
//cout << "else" << endl;
					list<int>::iterator it;
					for(it = map_reads.begin(); it != map_reads.end(); it++)
						data_tag[*it] = -1;
				}
*/
			} else { //right_id != -1
//cout << "11" << endl;
				bool fail_add = false;
				int length = forward_seq.length() - g_read_length - right_id.second;
				int start1 = i;
				string::size_type start2 = node_set[right_node].sequence.find(data[right_id.first]);
				if (start2 == string::npos)
					fail_add = true;
				string bubble_seq;
				if (length > 0)
					bubble_seq = forward_seq.substr(g_read_length, length);
				if ( p == right_node) {
					if ((!fail_add) && start2 <= start1)
						fail_add = true;
					int distance = 0;
					if (!fail_add)
						distance = start2 - start1 - g_read_length;
					//if ((!fail_add) && length == distance &&  distance > 0 && is_seq_similar(bubble_seq, node_set[p].sequence.substr(start2+g_read_length, distance)))
						//fail_add = true;
/*
					if ((!fail_add) && start1 < start2 && distance <= 0 && length > 0) {
						string sequence = node_set[p].sequence.substr(0, start1+g_read_length) + bubble_seq + node_set[p].sequence.substr(start2);
						node_set[p].sequence = sequence;
						set_reads_pos_in_node(kmer_hash, p);
						set_check_pos(p);
						forward_extend_by_coverage(kmer_hash, p);
						break_while = true;
						break;
					}
*/
					if (distance <= 0)
						fail_add = true;
					if (distance + length <= 4)
						fail_add = true;
					if (!fail_add) {
//cout << "p==q"<< endl;
//cout << i << endl;
						Node node1, node2, node3;
						node_id n1 = -1;
						node_id n2 = -1;
						node_id n3 = -1;
						if (length > 0) {
							node1.sequence = bubble_seq;
							n1 = add_node(node1);
							set_reads_pos_in_node(kmer_hash, n1);
							set_check_pos(n1);
						}
						node2.sequence = node_set[p].sequence.substr(0, start1+g_read_length);
						node3.sequence = node_set[p].sequence.substr(start1+g_read_length, distance);
						node_set[p].sequence = node_set[p].sequence.substr(start2);
						n2 = add_node(node2);
						n3 = add_node(node3);
						for (int j = 0; j < node_set[p].parents.size(); j++) {
							int p_i = node_set[p].parents[j];
							node_set[p_i].delete_child(p);
							node_set[p_i].add_child(n2);
						}
						node_set[n2].parents = node_set[p].parents;
						node_set[p].parents.clear();
						node_set[n2].add_child(n3);
						node_set[n3].add_parent(n2);
						node_set[n3].add_child(p);
						node_set[p].add_parent(n3);
						if (length > 0) {
							node_set[n1].add_child(p);
							node_set[p].add_parent(n1);
							node_set[n2].add_child(n1);
							node_set[n1].add_parent(n2);
						} else {
							node_set[n2].add_child(p);
							node_set[p].add_parent(n2);
						}
						set_reads_pos_in_node(kmer_hash, p);
						set_reads_pos_in_node(kmer_hash, n2);
						set_reads_pos_in_node(kmer_hash, n3);
						set_check_pos(p);
						set_check_pos(n2);
						set_check_pos(n3);
						forward_extend_by_coverage(kmer_hash, n3);
						forward_extend_by_coverage(kmer_hash, p);
						//if (length > 0)//
							//forward_extend_by_coverage(kmer_hash,n1);//
						break_while = true;
						break;
					}

				} else { //p != right_node
//cout<< "p1q" << endl;
//cout << i << endl;
					set<node_id> checked;
					if (right_node == 0 || is_circle(p, right_node, checked))
						fail_add = true;
							
					if (!fail_add) {
						Node node1, node2, node3;
						node_id n1 = -1;
						node_id n2 = -1;
						node_id n3 = -1;
						if (length > 0) {
							node1.sequence = bubble_seq;
							n1 = add_node(node1);
							set_reads_pos_in_node(kmer_hash, n1);
							set_check_pos(n1);
						}
						if (start1+g_read_length < node_set[p].sequence.length()) {
							node2.sequence = node_set[p].sequence.substr(0, start1+g_read_length);
							n2 = add_node(node2);
							set_reads_pos_in_node(kmer_hash, n2);
							set_check_pos(n2);
							node_set[p].sequence = node_set[p].sequence.substr(start1+g_read_length);
							set_reads_pos_in_node(kmer_hash, p);
							set_check_pos(p);
							for (int j = 0; j < node_set[p].parents.size(); j++) {
								int p_i = node_set[p].parents[j];
								node_set[p_i].delete_child(p);
								node_set[p_i].add_child(n2);
							}
							node_set[n2].parents = node_set[p].parents;
							node_set[p].parents.clear();
							node_set[n2].add_child(p);
							node_set[p].add_parent(n2);
							if (length > 0) {
								node_set[n2].add_child(p);
								node_set[p].add_parent(n2);
							}
						} else {
							if (length > 0) {
								node_set[p].add_child(n1);
								node_set[n1].add_parent(p);
							}
						}
						if (start2 > 0) {
							node3.sequence = node_set[right_node].sequence.substr(0, start2);
							n3 = add_node(node3);
							set_reads_pos_in_node(kmer_hash, n3);
							set_check_pos(n3);
							node_set[right_node].sequence = node_set[right_node].sequence.substr(start2);
							set_reads_pos_in_node(kmer_hash, right_node);
							set_check_pos(right_node);
							for (int j = 0; j < node_set[right_node].parents.size(); j++) {
								int p_i = node_set[right_node].parents[j];
								node_set[p_i].delete_child(right_node);
								node_set[p_i].add_child(n3);
							}
							node_set[n3].parents = node_set[right_node].parents;
							node_set[right_node].parents.clear();
							node_set[n3].add_child(right_node);
							node_set[right_node].add_parent(n3);									
						}
						if (length > 0) {
							node_set[right_node].add_parent(n1);
							node_set[n1].add_child(right_node);
						} else {
							if (start1+g_read_length < node_set[p].sequence.length()) {
								node_set[n2].add_child(right_node);
								node_set[right_node].add_parent(n2);
							} else {
								node_set[p].add_child(right_node);
								node_set[right_node].add_parent(p);
							}
						}
						forward_extend_by_coverage(kmer_hash, p);
						//if (length > 0)//
							//forward_extend_by_coverage(kmer_hash, n1);//
						break_while = true;
						break;
					}
				}


				if (fail_add) {
					if (forward_seq.length() > g_min_exon_length + g_read_length && pair_support(map_reads, p)){ //&& (!is_seq_similar(forward_seq.substr(g_read_length), node_set[p].sequence.substr(i + g_read_length))) ){


//cout << "failadd" << endl;
//cout << node_set[p].sequence.length() << endl;
//cout << forward_seq.length() << endl;

						
						Node node2 = node_set[p].sequence.substr(0, i + g_read_length);
						int n2 = add_node(node2);
						node_set[p].sequence = node_set[p].sequence.substr(i + g_read_length);
						for (int j = 0; j < node_set[p].parents.size(); j++) {
							int p_i = node_set[p].parents[j];
							node_set[p_i].delete_child(p);
							node_set[p_i].add_child(n2);
						}
						node_set[n2].parents = node_set[p].parents;
						node_set[p].parents.clear();
						node_set[p].add_parent(n2);
						node_set[n2].add_child(p);
						set_reads_pos_in_node(kmer_hash, n2);
						set_check_pos(n2);
						if (node_set[p].sequence.length() < g_min_exon_length && node_set[p].children.size() == 0) {

//cout << "if" << endl;
//cout << i << endl;
							set_reads_tag(kmer_hash, node_set[p].sequence, -1);
							node_set[p].sequence = forward_seq.substr(g_read_length);
							set_reads_pos_in_node(kmer_hash, p);
							set_check_pos(p);
							forward_extend_by_coverage(kmer_hash, p);
							break_while = true;
							break;
						} else {

//cout << "else" << endl;
//cout << i << endl;
							Node node1;
							node1.sequence = forward_seq.substr(g_read_length);
							int n1 = add_node(node1);
							node_set[n2].add_child(n1);
							node_set[n1].add_parent(n2);							
							set_reads_pos_in_node(kmer_hash, n1);
							set_reads_pos_in_node(kmer_hash, p);
							set_check_pos(p);
							set_check_pos(n1);
							forward_extend_by_coverage(kmer_hash, p);
							//forward_extend_by_coverage(kmer_hash, n1);//
							break_while = true;
							break;
						}
						
					}/* else {
						list<int>::iterator it;
						for(it = map_reads.begin(); it != map_reads.end(); it++)
							data_tag[*it] = -1;
					}
*/

				} //if (fail_add)

			}

		} // for

		if (break_while)
			break;

	}

}



void SplicingGraph::reverse_extend_by_coverage(KmerHash& kmer_hash, node_id p) {

	if (p < 0 || p >= node_sum)
		return;
cout << "reverse  check " << p << endl;
	while (node_set[p].reverse_check_pos.size() > 0) {

		bool break_while = false;
		pair<int, int> check_pos = node_set[p].reverse_check_pos.front();
		node_set[p].reverse_check_pos.pop_front();
cout << node_set[p].sequence.length() << " " << check_pos.first << " " << check_pos.second << endl;
		for (int i = check_pos.first; i < check_pos.second; i++) {
			string seed = node_set[p].sequence.substr(i, g_read_length);
			list<int> map_reads;
			string reverse_seq = reverse_extend(kmer_hash, seed, map_reads);
			list<int>::iterator it;
			for(it = map_reads.begin(); it != map_reads.end(); it++)
				data_tag[*it] = -1;
			if ((!pair_support(map_reads,p)) || reverse_seq.length() < g_read_length + g_min_exon_length) {
				/*list<int>::iterator it;
				for(it = map_reads.begin(); it != map_reads.end(); it++)
					data_tag[*it] = -1;
*/
				continue;
			}
			if (!is_seq_similar(reverse_seq, node_set[p].sequence.substr(0, i + g_read_length), 'R')) {
//cout << i << endl;
				Node node1, node2;
				node1.sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
				node_id n1 = add_node(node1);
				set_reads_pos_in_node(kmer_hash, n1);
				set_check_pos(n1);
				node2.sequence = node_set[p].sequence.substr(0, i);
				node_id n2 = add_node(node2);
				node_set[p].sequence = node_set[p].sequence.substr(i);
				set_reads_pos_in_node(kmer_hash, n2);
				set_reads_pos_in_node(kmer_hash, p);
				set_check_pos(n2);
				set_check_pos(p);
				for (int j = 0; j < node_set[p].parents.size(); j++) {
					int p_i = node_set[p].parents[j];
					node_set[p_i].delete_child(p);
					node_set[p_i].add_child(n2);
				}
				node_set[n2].parents = node_set[p].parents;
				node_set[p].parents.clear();
				node_set[n2].add_child(p);
				node_set[p].add_parent(n2);
				node_set[n1].add_child(p);
				node_set[p].add_parent(n1);
				reverse_extend_by_coverage(kmer_hash, p);
				//reverse_extend_by_coverage(kmer_hash, n1);//
				break_while = true;
				break;

			}/* else {
				list<int>::iterator it;
				for(it = map_reads.begin(); it != map_reads.end(); it++)
					data_tag[*it] = -1;
			}
*/

		}

		if (break_while)
			break;
	}

}


void SplicingGraph::set_check_pos(node_id p) {

//cout << "set check" << endl;
	node_set[p].forward_check_pos.clear();
	node_set[p].reverse_check_pos.clear();
	if (node_set[p].sequence.length() <= g_read_length + 1)
		return;
	if (node_set[p].cov_reads.size() < 10)
		return;

	int node_len = node_set[p].sequence.length();
	vector<float> node_cov(node_len, 0.0);
	for (int i = 0; i < node_set[p].cov_reads.size(); i++)
		for (int j = node_set[p].cov_reads[i].second; j < node_set[p].cov_reads[i].second+g_read_length; j++)
			node_cov[j] = node_cov[j] + 1.0;


	bool first_reverse = true;
	bool first_forward = true;
	int reverse_start = 0;
	int reverse_end = 0;
	int forward_start = 0;
	int forward_end = 0;
	bool forward_tag = false;
	bool reverse_tag = false;
	vector<pair<int, int> > candi_forward;
	vector<pair<int, int> > candi_reverse;
//for (int i = 0; i < node_cov.size(); i++)
//cout << node_cov[i] << ",";
//cout << endl;

	for (int i = 1; i < (int)node_cov.size()-g_read_length - 1; i++) {
//cout << i << endl;
		if (node_cov[i] / node_cov[i+g_min_same_len-1] < 0.5 || node_cov[i+g_min_same_len-1] - node_cov[i] > 100){

			if (forward_tag) {
				pair<int, int> fp;
				fp.first = forward_start;
				fp.second = forward_end;
				node_set[p].forward_check_pos.push_back(fp);
				first_forward = true;
				forward_tag = false;
			}

			if (first_reverse) {
				reverse_start = i;
				first_reverse = false;
			}
			reverse_end = i + g_min_same_len-1;
			reverse_tag = true;
		
		} else {

			if (node_cov[i+g_min_same_len-1] / node_cov[i] < 0.5 || node_cov[i] - node_cov[i+g_min_same_len-1] > 100) {

				if (reverse_tag) {
					pair<int,int> rp;
					rp.first = reverse_start;
					rp.second = reverse_end;
					node_set[p].reverse_check_pos.push_back(rp);
					first_reverse = true;
					reverse_tag = false;
				}

				if (first_forward) {
					forward_start = i;
					first_forward = false;
				}
				forward_end = i + g_min_same_len - 1;
				forward_tag = true;
			}else {

				if (forward_tag) {
					pair<int, int> fp;
					fp.first = forward_start;
					fp.second = forward_end;
					node_set[p].forward_check_pos.push_back(fp);
					first_forward = true;
					forward_tag = false;
				}

				if (reverse_tag) {
					pair<int,int> rp;
					rp.first = reverse_start;
					rp.second = reverse_end;
					node_set[p].reverse_check_pos.push_back(rp);
					first_reverse = true;
					reverse_tag = false;
				}				
			}

		}

	}
	
	if (candi_forward.size() > 0) {
		int n = 1;
		pair <int, int> current_pos = candi_forward[0];
		
		while (n < candi_forward.size()) {
			if (candi_forward[n].first <= current_pos.second) {
				current_pos.second = candi_forward[n].second;
			} else {
				node_set[p].forward_check_pos.push_back(current_pos);
				current_pos = candi_forward[n];
			}		
			n = n + 1;
		}

		node_set[p].forward_check_pos.push_back(current_pos);
	}

	if (candi_reverse.size() > 0) {
		int n = 1;
		pair<int, int> current_pos = candi_reverse[0];
		while(n < candi_reverse.size()) {
			if (candi_reverse[n].first < current_pos.second)
				current_pos.second = candi_reverse[n].second;
			else {
				node_set[p].reverse_check_pos.push_back(current_pos);
				current_pos = candi_reverse[n];
			}
			n = n + 1;
		}
		node_set[p].reverse_check_pos.push_back(current_pos);
	}
/*
pair<int, int> xxx;
xxx.first = 0;
xxx.second = node_set[p].sequence.length()-g_read_length;
node_set[p].reverse_check_pos.push_back(xxx);
node_set[p].forward_check_pos.push_back(xxx);
*/

}
