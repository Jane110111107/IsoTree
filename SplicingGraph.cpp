
#include "SplicingGraph.h"
#include<iostream>

using namespace std;


int get_compatible_len(const string& str1, const string& str2) {

	int compatible_len = 0;
	for (int i = g_min_same_len-1; i > 5; --i) {
		const string& comp1 = str1.substr(0, i);
		const string& comp2 = str2.substr(str2.length()-i);
		if (comp1 == comp2) {
			compatible_len = i;
			break;
		}
	}
	
	return compatible_len;

}


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

		for (int i = 0; i < length; ++i) {

			if (str1[i] != str2[i])
				mismatch++;

		}

	} else {

		for (int i = 0; i < length; ++i) {

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
		for(int i = 0; i < length; ++i) {	
			if(str1[i] != str2[i])						
				mismatch++;
		}
	} else {
		for (int i = 0; i < length; ++i) {
			if (str1[str1.length()-i-1] != str2[str2.length()-i-1])
				mismatch++;
		}
	}
	return (mismatch <= tag);

}


bool compatible(const string& str1, const string& str2) {

	for(unsigned int i = 0; i <= str2.length()-g_kmer_length; ++i) {

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

	for (int i = 0; i < node_sum; ++i) {

		if (!node_set[i].parents.empty())
			node_set[i].parents.clear();

	}

	for (size_t i = 0; i < node_sum; ++i) {

		vector<node_id>::iterator it;

		for (it = node_set[i].children.begin(); it != node_set[i].children.end(); ++it)
			node_set[*it].add_parent(i);

	}

}



int SplicingGraph::add_node(Node& node) {

	node_set.push_back(node);
	return (node_sum++);

}


void SplicingGraph::refine_trunk(ReadHash& read_hash, KmerHash& kmer_hash) {

//refine reverse
	int check_len;
	if (node_set[0].sequence.length() < (g_read_length + 2*g_kmer_length)*2)
		check_len = node_set[0].sequence.length() / 2;
	else
		check_len = 2*g_kmer_length;
	if (check_len + g_read_length > node_set[0].sequence.length())
		check_len = node_set[0].sequence.length() - g_read_length;
	if (check_len < 1)
		return;

	read_int_type read_int;
	for (int i = 0; i < check_len; ++i) {
		const string& seed = node_set[0].sequence.substr(i, g_read_length);
		read_int = get_read_int_type(seed);
		if (read_hash[read_int].size() < 1)
			continue;
		const string& sequence = reverse_extend(read_hash, kmer_hash, seed);
		if ((int)sequence.length() > g_read_length + g_min_exon_length && pair_support(read_hash, sequence, 0)) {
			set_reads_tag(read_hash, node_set[0].sequence.substr(0, i+g_read_length-1), -4);
			node_set[0].sequence = sequence + node_set[0].sequence.substr(i+g_read_length);
			set_reads_tag(read_hash, sequence, 0);
			break;
		} else { //recover reads
			set_reads_tag(read_hash, sequence, -4);
		}
	}

	if (node_set[0].sequence.length() < (g_read_length + 2*g_kmer_length)*2)
		check_len = node_set[0].sequence.length() / 2;
	else
		check_len = 2*g_kmer_length;
	if (check_len + g_read_length > node_set[0].sequence.length())
		check_len = node_set[0].sequence.length() - g_read_length;

	for (int i = node_set[0].sequence.length()-g_read_length-1; i >= node_set[0].sequence.length()-g_read_length-check_len; --i) {
		const string& seed = node_set[0].sequence.substr(i, g_read_length);
		const string& sequence = forward_extend(read_hash, kmer_hash, seed);
		if ((int)sequence.length() > g_read_length + g_min_exon_length && pair_support(read_hash, sequence, 0)) {
			set_reads_tag(read_hash, node_set[0].sequence.substr(i), -4);
			node_set[0].sequence = node_set[0].sequence.substr(0,i) + sequence;
			set_reads_tag(read_hash, sequence, 0);
			break;
		} else { //recover reads
			set_reads_tag(read_hash, sequence, -4);
		}
	}

}


bool SplicingGraph::is_trunk(ReadHash& read_hash, const string& trunk) {

	vector<int> read_set;
	int used_sum = 0;
	int total_sum = 0;
	for (int i = 0; i <= trunk.length()-g_read_length; ++i) {
		const string& read = trunk.substr(i, g_read_length);
		read_int_type read_int = get_read_int_type(read);
		read_set = read_hash[read_int];
		for (int j = 0; j < read_set.size(); ++j) {
			if (data_used[read_set[j]]) {
				++used_sum;
				break;
			}
		}
		if (read_set.size() > 0)
			++total_sum;
	}
	if (used_sum > 0.8*total_sum && trunk.length() < 1000)
		return false;
	if (used_sum > 0.9*total_sum)
		return false;

	return true;

/*
	vector<pair<int, int> > reads_pos;
	//reads_pos.reserve(50000);
	map_reads_to_sequence(read_hash, trunk, reads_pos);
	vector<pair<int, int> > pair_reads = reads_pos;
	pair_reads.clear();
	int mate_id = -1;
	for (int i =0; i < reads_pos.size(); ++i) {
		if (reads_pos[i].first >= max_read_id && data_tag[reads_pos[i].first-max_read_id] > -2)
			pair_reads.push_back(reads_pos[i]);
		if (reads_pos[i].first < max_read_id && data_tag[reads_pos[i].first+max_read_id] > -2)
			pair_reads.push_back(reads_pos[i]);
	}
	if (pair_reads.size() < 5) {
		return false;
	}
	int total_cov_len = 0;
	int current_cov_len = g_read_length;
	int current_cov_start = pair_reads[0].second;
	bool add_tag = true;
	for (int i = 0; i < pair_reads.size()-1; ++i) {
		if (pair_reads[i].second + g_read_length >= pair_reads[i+1].second) {
			add_tag = true;
			current_cov_len = pair_reads[i+1].second + g_read_length - current_cov_start;
		} else {
			add_tag = false;
			total_cov_len = total_cov_len + current_cov_len;
			current_cov_len = g_read_length;
			current_cov_start = pair_reads[i+1].second;
		}
	}
	if (add_tag = true)
		total_cov_len = total_cov_len + current_cov_len;

	if (total_cov_len < trunk.length()*0.6)
		return false;
	return true;
*/
	
}

void SplicingGraph::set_trunk_check_pos(ReadHash& read_hash, node_id p) {

	int node_len = node_set[p].sequence.length();
	vector<float> node_cov(node_len, 0.0);
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i)
		for (int j = node_set[p].cov_reads[i].second; j < node_set[p].cov_reads[i].second+g_read_length; ++j)
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
	for (int i = g_read_length; i < node_cov.size()-2*g_read_length; ++i) {

		if (node_cov[i] < node_cov[i+g_min_same_len-1] * 0.4){

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

			if (node_cov[i+g_min_same_len-1] < node_cov[i] * 0.4) {

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
}



void SplicingGraph::refine_forward_trunk(ReadHash& read_hash, KmerHash& kmer_hash) {

	if (node_set[0].sequence.length() <= g_read_length)
		return;
	const string& trunk = node_set[0].sequence;
	list<int> map_reads;
	int trunk_len = trunk.length();
	int check_len = g_min_same_len*2;
	if (check_len + g_read_length > trunk_len)
		check_len = trunk_len - g_read_length;

	for (int i = 1; i <= check_len; ++i) {
		const string& read = trunk.substr(trunk_len-i-g_read_length, g_read_length);
		read_int_type read_int = get_read_int_type(read);
		if (read_hash[read_int].size() == 0)
			continue;
		const string& forward_trunk = forward_extend(read_hash, kmer_hash, read);
		if (forward_trunk.length() >= g_read_length + g_min_exon_length) {
			node_set[0].sequence = trunk.substr(0, trunk_len-i-g_read_length) + forward_trunk;
			if (trunk_len-i-2*g_read_length+1 <= 0)
				set_reads_tag(read_hash, trunk, -4);
			else {
				const string& recover_str = trunk.substr(trunk_len-i-2*g_read_length+1);
				set_reads_tag(read_hash, recover_str, -4);
			}
			break;
		} else {
			set_reads_tag(read_hash,forward_trunk,-4);
		}
	}

}


void SplicingGraph::refine_reverse_trunk(ReadHash& read_hash, KmerHash& kmer_hash) {

	if (node_set[0].sequence.length() <= g_read_length)
		return;
	const string& trunk = node_set[0].sequence;
	int trunk_len = trunk.length();
	int check_len = g_min_same_len*2;
	if (check_len + g_read_length > trunk_len)
		check_len = trunk_len - g_read_length;
		
	for (int i = 1; i <= check_len; ++i) {

		const string& read = trunk.substr(i, g_read_length);
		read_int_type read_int = get_read_int_type(read);
		if (read_hash[read_int].size() == 0)
			continue;
		const string& reverse_trunk = reverse_extend(read_hash, kmer_hash, read);
		if (reverse_trunk.length() >= g_read_length + g_min_exon_length) {
			node_set[0].sequence = reverse_trunk + trunk.substr(i+g_read_length);
			if (i + 2*g_read_length-1 > trunk_len)
				set_reads_tag(read_hash, trunk, -4);
			else {
				const string& recover_str = trunk.substr(0, i+g_read_length);
				set_reads_tag(read_hash, recover_str, -4);
			}
			break;
		} else {
			set_reads_tag(read_hash,reverse_trunk,-4);
		}
	}

}




bool SplicingGraph::refine_reverse_by_pair(ReadHash& read_hash, KmerHash& kmer_hash, node_id p) {

//cout << "refine reverse by pair ..." << endl;
	if (node_set[p].sequence.length() < g_read_length+g_refine_check_len)
		return false;
	bool return_tag = false;
	vector<int> check_reads;
	int extend_time = 0;
	while (1) {		
		const string& check_seq = node_set[p].sequence.substr(0, g_read_length+g_refine_check_len);
		check_reads.clear();
		set_reads_tag(read_hash, check_seq,check_reads, p);
		if (check_reads.size() == 0)
			break;
		bool is_add = false;
		for (int i = 0; i < check_reads.size(); ++i){
			int read_id = check_reads[i];
			if (read_id >= max_read_id && data_tag[read_id-max_read_id] < -2){
				++extend_time;
				if (extend_time > 20)
					break;
				string extend_seq = data[read_id-max_read_id];
				data_tag[read_id-max_read_id] = -1;
				const string& stop_seq = node_set[p].sequence.substr(0, g_min_same_len-1);
				is_add = forward_extend(read_hash, kmer_hash,extend_seq,stop_seq);
				if (is_add) {
					return_tag = true;
					string seq = reverse_extend(read_hash, kmer_hash, extend_seq);				
					node_set[p].sequence = seq + node_set[p].sequence;
					break;
				} else {
					if (extend_seq.length() > 1000) {
						Node node1;
						node1.sequence = extend_seq;
						node_id n1 = add_node(node1);
						set_reads_tag(read_hash, extend_seq, n1);
					} else {
						set_reads_tag(read_hash,extend_seq,-4);
					}
				}
			}

		}
		if (!is_add)
			break;
	}

	return return_tag;

}





bool SplicingGraph::refine_forward_by_pair(ReadHash& read_hash, KmerHash& kmer_hash, node_id p) {

//cout << "refine forward by pair..." << endl;
	if (node_set[p].sequence.length() < g_read_length+g_refine_check_len)
		return false;
	bool return_tag = false;
	vector<int> check_reads;
	int extend_time = 0;
	while (1) {		
		const string& check_seq = node_set[p].sequence.substr(node_set[p].sequence.length()-g_read_length-g_refine_check_len);
		check_reads.clear();
		set_reads_tag(read_hash, check_seq,check_reads, p);
		if (check_reads.size() == 0)
			break;
		bool is_add = false;
		for (int i = 0; i < check_reads.size(); ++i){
			int read_id = check_reads[i];
			if (read_id < max_read_id && data_tag[read_id+max_read_id] < -2){
				++extend_time;
				if (extend_time > 20)
					break;
				data_tag[read_id+max_read_id] = -1;
				string extend_seq = data[read_id+max_read_id];
				const string& stop_seq = node_set[p].sequence.substr(node_set[p].sequence.length()-g_min_same_len+1);
				is_add = reverse_extend(read_hash, kmer_hash,extend_seq, stop_seq);
				if (is_add) {
					return_tag = true;
					string seq = forward_extend(read_hash, kmer_hash, extend_seq);				
					node_set[p].sequence = node_set[p].sequence + seq;
					break;
				} else {
					if (extend_seq.length() > 1000) {
						Node node1;
						node1.sequence = extend_seq;
						node_id n1 = add_node(node1);
						set_reads_tag(read_hash, extend_seq, n1);
					} else {
						set_reads_tag(read_hash,extend_seq,-4);
					}
				}
			}

		}
		if (!is_add)
			break;
	}

	return return_tag;

}





bool SplicingGraph::build(ReadHash& read_hash, KmerHash& kmer_hash, int seed, vector<pair<string,float> >& transcripts) {

	data_tag[seed] = -1;
	const string& seed_str = data[seed];

	read_int_type read_int = get_read_int_type(seed_str);
	vector<int> read_set = read_hash[read_int];
	for (int i = 0; i < read_set.size(); ++i)
		data_tag[read_set[i]] = -1;
	
	cout << " Get trunk..." << endl;
	const string& left = reverse_extend(read_hash, kmer_hash, seed_str);
	const string& right = forward_extend(read_hash, kmer_hash, seed_str);
	string trunk = left;
	if (right.length() > g_read_length)
		trunk = left + right.substr(g_read_length);

	if(static_cast<int>(trunk.length()) < g_min_trunk_length) {
		cout << "The trunk is too short!" << endl;
		set_reads_tag(read_hash, trunk, -4);		
		return false;
	}
/*
	if (!is_trunk(read_hash, trunk)) {
		set_reads_tag(read_hash, trunk, -4);		
		return false;
	}
*/


	Node node(trunk);
	int p = add_node(node);
	//refine_forward_trunk(read_hash, kmer_hash);
	//refine_reverse_trunk(read_hash, kmer_hash);
	set_reads_pos_in_node(read_hash, p);
//*
	if (g_is_paired_end) {
		if (refine_reverse_by_pair(read_hash,kmer_hash,p) || refine_forward_by_pair(read_hash,kmer_hash,p))
			set_reads_pos_in_node(read_hash, p);
	}
//*/
	if (g_mode == 2) {
		//set_trunk_check_pos(read_hash, p);
		cout << "branch check and extend..." << endl;
		branch_extend_by_coverage(read_hash, kmer_hash);
		cout << "refine graph..." << endl;
		if (g_is_paired_end)
			refine_graph(read_hash, kmer_hash);
		//trim_graph(read_hash);
	}
		topological_sort();
		get_transcripts(read_hash, transcripts);

	if (transcripts.size() == 0)
		return false;

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



void SplicingGraph::refine_graph(ReadHash& read_hash, KmerHash& kmer_hash) {

	int num = 0;
	for (int i = 0; i < node_sum; ++i) {
		bool is_forward = false;
		bool is_reverse = false;
		if (node_set[i].children.size() == 0) {
			num = 0;
			is_forward=refine_forward_by_pair(read_hash, kmer_hash, i);
		}
		if (node_set[i].parents.size() == 0) {
			num = 0;
			is_reverse=refine_reverse_by_pair(read_hash, kmer_hash, i);
		}
		if (is_forward || is_reverse)
			set_reads_pos_in_node(read_hash, i);

	}

}




void SplicingGraph::refine_forward(ReadHash& read_hash, KmerHash& kmer_hash, node_id p, int& num) {

	if (p < 0 || p >= node_sum)
		return;
	if (node_set[p].sequence.length() <= g_read_length)
		return;
//cout << "refine forward..." << endl;

	bool is_add = false;
	int check_len = g_kmer_length;
	if (node_set[p].sequence.length() < check_len + g_read_length)
		check_len = node_set[p].sequence.length()-g_read_length;

	vector<int> read_set;
	read_int_type read_int;
	int mate_id;
	int check_times = 0;
	const string& stop_seq = node_set[p].sequence.substr(node_set[p].sequence.length()-g_kmer_length+1);
	for (int i = node_set[p].sequence.length()-g_read_length-1; i > node_set[p].sequence.length()-g_read_length-check_len; --i) {//forward check
		const string& seed_read = node_set[p].sequence.substr(i, g_read_length);
		read_int = get_read_int_type(seed_read);
		read_set = read_hash[read_int];
		for (int j = 0; j < read_set.size(); ++j) {
			if (read_set[j] < max_read_id && data_tag[read_set[j]+max_read_id] <-3) {
				++check_times;
				string mate_read = data[read_set[j]+max_read_id];
				data_tag[read_set[j]+max_read_id] = -1;
				is_add = reverse_extend(read_hash, kmer_hash,mate_read, stop_seq) ;
				string::size_type start = mate_read.find(stop_seq);
				if (is_add) {
					const string& extend_seq = forward_extend(read_hash, kmer_hash, mate_read);
					node_set[p].sequence = node_set[p].sequence + extend_seq.substr(start+g_kmer_length);
					set_reads_tag(read_hash, extend_seq, p);
					++num;
					if (num > 2)
						return;
					else {
						refine_forward(read_hash, kmer_hash, p, num);
						return;
					}
				} else {
					string check_seq = stop_seq.substr(1);
					while (check_seq.length() > 10) {
						start = mate_read.find(check_seq);
						if (start != string::npos) {
							const string& extend_seq = forward_extend(read_hash, kmer_hash, mate_read);
							node_set[p].sequence = node_set[p].sequence + extend_seq.substr(start+check_seq.length());
							set_reads_tag(read_hash, extend_seq, p);
							++num;
							if (num > 2)
								return;
							else {
								refine_forward(read_hash, kmer_hash, p, num);
								return;
							}
						}
						check_seq = check_seq.substr(1);
					}
				}
				set_reads_tag(read_hash, mate_read, -4);
				if (check_times > 50)
					return;
			}

		}
		
	}

	return;

}



void SplicingGraph::refine_reverse(ReadHash& read_hash, KmerHash& kmer_hash, node_id p, int& num) {

	if (p < 0 || p >= node_sum)
		return;
	if (node_set[p].sequence.length() <= g_read_length)
		return;
//cout << "refine reverse..." << endl;
	//set_reads_pos_in_node(kmer_hash, p);
	bool is_add = false;
	int check_times = 0;
	int check_len = g_kmer_length;
	if (node_set[p].sequence.length() < check_len + g_read_length)
		check_len = node_set[p].sequence.length()-g_read_length;
	const string& stop_seq = node_set[p].sequence.substr(0, g_kmer_length-1);
	vector<int> read_set;
	read_int_type read_int;
	int mate_id;

	for (int i = 0; i < check_len; ++i){//reverse check
		const string& seed_read = node_set[p].sequence.substr(i, g_read_length);
		read_int = get_read_int_type(seed_read);
		read_set = read_hash[read_int];
		for (int j = 0; j < read_set.size(); ++j) {
			if (read_set[j] >= max_read_id && data_tag[read_set[j]-max_read_id] < -3){
				++check_times;
				string mate_read = data[read_set[j]-max_read_id];
				data_tag[read_set[j]-max_read_id] = -1;
				is_add = forward_extend(read_hash, kmer_hash,mate_read, stop_seq) ;
				string::size_type start = mate_read.find(stop_seq);
				if (is_add) {
					int same_len = mate_read.length() - start;
					const string& extend_seq = reverse_extend(read_hash, kmer_hash, mate_read);
					node_set[p].sequence =  extend_seq.substr(extend_seq.length()-same_len) + node_set[p].sequence;
					set_reads_tag(read_hash, extend_seq, p);
					++num;
					if (num > 2)
						return;
					else {
						refine_reverse(read_hash, kmer_hash, p, num);
						return;
					}
				} else {
					string check_seq = stop_seq.substr(0, stop_seq.length()-1);
					while (check_seq.length() > 10) {
						start = mate_read.find(check_seq);
						if (start != string::npos) {
							int same_len = mate_read.length() - start;
							const string& extend_seq = reverse_extend(read_hash, kmer_hash, mate_read);
							node_set[p].sequence =  extend_seq.substr(extend_seq.length()-same_len) + node_set[p].sequence;
							set_reads_tag(read_hash, extend_seq, p);
							++num;
							if (num > 2)
								return;
							else {
								refine_reverse(read_hash, kmer_hash, p, num);	
								return;
							}
						}
						check_seq = check_seq.substr(0, check_seq.length()-1);
					}
				}
				set_reads_tag(read_hash, mate_read, -4);
				if (check_times > 50)
					return; 
			}
		}
		
	}

	return;
}


bool SplicingGraph::pair_support(ReadHash& read_hash, const string& sequence, node_id p) {

//cout << "pair support " << endl;
	if (p < 0 || p >= node_sum)
		return false;

	if (!g_is_paired_end)
		return true;

	int support = 0;
	vector<int> read_set;
	read_int_type read_int;
	int mate_id;
	for (int i = 0; i < sequence.length()-g_read_length; ++i) {
		const string& read_seq = sequence.substr(i, g_read_length);
		read_int = get_read_int_type(read_seq);
		read_set = read_hash[read_int];
		for (int j = 0; j < read_set.size(); ++j) {
			if (read_set[j] >= max_read_id)
				mate_id = read_set[j] - max_read_id;
			else
				mate_id = read_set[j] + max_read_id;
			if (data_tag[mate_id] == p) {
				++support;
				if (support >= g_min_reads_support_branch)
					return true;
			}
		}
	}
	
	return false;

}


string SplicingGraph::forward_extend(ReadHash& read_hash, KmerHash& kmer_hash, const string& seed_contig, int color_tag) {

//cout << "forward extend..." << endl;
	if (seed_contig.length() < g_read_length)
		return seed_contig;

	string contig = seed_contig;
//*
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	int mate_id, r_d, candi_id;
	kmer_int_type tag_kmer, kmer_int;
	vector<int> seed_readset;
	vector<int> tag_readset;
	while (1) {

		const string& forward_kmer = contig.substr(contig.length()-g_kmer_length, g_kmer_length);
		kmer_int = kmer_to_int(forward_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(contig.length() - g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		candi_id = -1;
		string::size_type start;
		int s_c = 0;
		if (g_is_paired_end) {
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				if (r_d >= max_read_id)
					mate_id = r_d - max_read_id;
				else
					mate_id = r_d + max_read_id;
				if (data_tag[mate_id] > -2) {
					start = data[r_d].find(forward_kmer);
					if (start != string::npos && start >= s_c && data[r_d].substr(0,start+g_kmer_length) == contig.substr(contig.length()-start-g_kmer_length)) {
						s_c = start;
						candi_id = r_d;
						if (start == g_read_length-g_kmer_length-1)
							break;
					}
				}
			}
		}
		if (candi_id > -1) {
			contig = contig + data[candi_id].substr(s_c + g_kmer_length);
			read_int = get_read_int_type(data[candi_id]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
			continue;
		}
		}
//*
		s_c = 0;
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				start = data[r_d].find(forward_kmer);
				if (start != string::npos && data[r_d].substr(0,start+g_kmer_length) == contig.substr(contig.length()-start-g_kmer_length)) {
					read_int = get_read_int_type(data[r_d]);
					if (start > s_c) {
						s_c = start;
						candi_id = r_d;
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
					} else {
						if (start == s_c && abs((int)read_hash[read_int].size() - ref_cov) < current_cov) {	
							s_c = start;
							candi_id = r_d;
							current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						}
					}
				}
			}
		}

		if (candi_id > -1) {
			contig = contig + data[candi_id].substr(s_c + g_kmer_length);
			read_int = get_read_int_type(data[candi_id]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
		} else {
			break;
		}
	}
/*
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	int mate_id, r_d, candi_id, candi_id1, s_c, s_c1;
	kmer_int_type tag_kmer, kmer_int;
	vector<int> seed_readset;
	vector<int> tag_readset;
	while (1) {

		const string& forward_kmer = contig.substr(contig.length()-g_kmer_length, g_kmer_length);
		kmer_int = kmer_to_int(forward_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(contig.length() - g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		candi_id = -1;
		candi_id1 = -1;
		string::size_type start;
		s_c = 0;
		s_c1 = 0;
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				start = data[r_d].find(forward_kmer);
				if (start == string::npos)
					continue;
				if (data[r_d].substr(0,start+g_kmer_length) == contig.substr(contig.length()-start-g_kmer_length)) {
					if (r_d >= max_read_id && data_tag[r_d-max_read_id] > -2 && start >= s_c1) {
						s_c1 = start;
						candi_id1 = r_d;
						if (start == g_read_length-g_kmer_length-1)
							break;
					}
					read_int = get_read_int_type(data[r_d]);
					if (start > s_c) {
						s_c = start;
						candi_id = r_d;
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
					} else {
						if (start == s_c && abs((int)read_hash[read_int].size() - ref_cov) < current_cov) {	
							s_c = start;
							candi_id = r_d;
							current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						}
					}
				}
			}
		}

		if (candi_id1 > -1) {
			contig = contig + data[candi_id1].substr(s_c1 + g_kmer_length);
			read_int = get_read_int_type(data[candi_id1]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
		} else {
			if (candi_id > -1) {
				contig = contig + data[candi_id].substr(s_c + g_kmer_length);
				read_int = get_read_int_type(data[candi_id]);
				candi_set = read_hash[read_int];
				for (int j = 0; j < candi_set.size(); ++j) 
					data_tag[candi_set[j]] = color_tag;
			} else {
				break;
			}
		}
	}
*/

/*
	string find_seq;
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	int mate_id, r_d, candi_id;
	kmer_int_type tag_kmer, kmer_int;
	vector<int> seed_readset;
	vector<int> tag_readset;
	while (1) {

		const string& forward_kmer = contig.substr(contig.length()-g_kmer_length, g_kmer_length);
		kmer_int = kmer_to_int(forward_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(contig.length() - g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		string candi;
		candi_id = -1;

		for (int i = 0; i < seed_readset.size(); ++i)
			data_find[seed_readset[i]] = true;

		for (int i = g_max_same_len; i >= g_min_same_len; --i) {
			const string& first_kmer = contig.substr(contig.length()-i,g_kmer_length);
			tag_kmer = kmer_to_int(first_kmer);	 
	       	 	tag_readset = kmer_hash[tag_kmer];
			const string& ref_seq = contig.substr(contig.length()-i);
			vector<int>::iterator b = tag_readset.begin();
	        	while (b != tag_readset.end()) {
				find_seq = data[*b].substr(0,i);
				if (data_find[*b] && data_tag[*b] < -2 && find_seq == ref_seq) {
					if (*b >= max_read_id)
						mate_id = *b - max_read_id;
					else
						mate_id = *b + max_read_id;
					if (data_tag[mate_id] > -2) {
						candi = data[*b].substr(i);
						candi_id = *b;
						break;
					}
					read_int = get_read_int_type(data[*b]);
					if (abs((int)read_hash[read_int].size() - ref_cov) < current_cov) {
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						candi = data[*b].substr(i);
						candi_id = *b;
					}
				}
				++b;				
			}//while(*b)
			if (candi_id > -1) {
				contig = contig + candi;
				read_int = get_read_int_type(data[candi_id]);
				candi_set = read_hash[read_int];
				for (int j = 0; j < candi_set.size(); ++j) {
					data_tag[candi_set[j]] = color_tag;
				}
				break;
			} 		
		}//for i
		if (candi_id == -1)
			break;
	} //while
*/
	return contig;
}



bool SplicingGraph::forward_extend(ReadHash& read_hash, KmerHash& kmer_hash, string& contig, const string& stop_seq, int color_tag) {

//cout << "forward extend..." << endl;
	if (contig.length() < g_read_length)
		return false;


//*
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	int mate_id, r_d, candi_id;
	kmer_int_type tag_kmer, kmer_int;
	vector<int> seed_readset;
	vector<int> tag_readset;
	while (1) {

		const string& forward_kmer = contig.substr(contig.length()-g_kmer_length, g_kmer_length);
		kmer_int = kmer_to_int(forward_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(contig.length() - g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		candi_id = -1;
		string::size_type start;
		int s_c = 0;
		if (g_is_paired_end) {
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				if (r_d >= max_read_id)
					mate_id = r_d - max_read_id;
				else
					mate_id = r_d + max_read_id;
				if (data_tag[mate_id] > -2) {
					start = data[r_d].find(forward_kmer);
					if (start != string::npos && start >= s_c && data[r_d].substr(0,start+g_kmer_length) == contig.substr(contig.length()-start-g_kmer_length)) {
						s_c = start;
						candi_id = r_d;
						if (start == g_read_length-g_kmer_length-1)
							break;
						int compatible_len = get_compatible_len(stop_seq, data[r_d].substr(g_read_length-g_min_same_len+1));
						if (compatible_len != 0){
							contig = contig.substr(0, contig.length()-i) + data[r_d].substr(0, g_read_length-compatible_len);
							data_tag[r_d] = color_tag;
							return true;
						}
					}
				}
			}
		}
		if (candi_id > -1) {
			contig = contig + data[candi_id].substr(s_c + g_kmer_length);
			read_int = get_read_int_type(data[candi_id]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
			continue;
		}
		}
//*
		s_c = 0;
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				start = data[r_d].find(forward_kmer);
				if (start != string::npos && data[r_d].substr(0,start+g_kmer_length) == contig.substr(contig.length()-start-g_kmer_length)) {
					int compatible_len = get_compatible_len(stop_seq, data[r_d].substr(g_read_length-g_min_same_len+1));
					if (compatible_len != 0){
						contig = contig.substr(0, contig.length()-i) + data[r_d].substr(0, g_read_length-compatible_len);
						data_tag[r_d] = color_tag;
						return true;
					}
					read_int = get_read_int_type(data[r_d]);
					if (start > s_c) {
						s_c = start;
						candi_id = r_d;
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
					} else {
						if (start == s_c && abs((int)read_hash[read_int].size() - ref_cov) < current_cov) {	
							s_c = start;
							candi_id = r_d;
							current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						}
					}
				}
			}
		}

		if (candi_id > -1) {
			contig = contig + data[candi_id].substr(s_c + g_kmer_length);
			read_int = get_read_int_type(data[candi_id]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
		} else {
			break;
		}
	}


/*
	string find_seq;
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	int mate_id, r_d, candi_id;
	kmer_int_type tag_kmer, kmer_int;
	vector<int> seed_readset;
	vector<int> tag_readset;
	string::size_type start, start1;
	while (1) {
		const string& forward_kmer = contig.substr(contig.length()-g_kmer_length, g_kmer_length);
		kmer_int = kmer_to_int(forward_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(contig.length() - g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		string candi;
		candi_id = -1;
		for (int i = 0; i < seed_readset.size(); ++i)
			data_find[seed_readset[i]] = true;
		for (int i = g_max_same_len; i >= g_min_same_len; --i) {
			const string& first_kmer = contig.substr(contig.length()-i,g_kmer_length);
			tag_kmer = kmer_to_int(first_kmer);	 
	       	 	tag_readset = kmer_hash[tag_kmer];
			const string& ref_seq = contig.substr(contig.length()-i);
			vector<int>::iterator b = tag_readset.begin();
	        	while (b != tag_readset.end()) {
				find_seq = data[*b].substr(0,i);
				if (data_find[*b] && data_tag[*b] < -2 && find_seq == ref_seq) {

					int compatible_len = get_compatible_len(stop_seq, data[*b].substr(g_read_length-g_min_same_len+1));
					if (compatible_len != 0){
						contig = contig.substr(0, contig.length()-i) + data[*b].substr(0, g_read_length-compatible_len);
						data_tag[*b] = color_tag;
						return true;
					}
					if (*b >= max_read_id)
						mate_id = *b - max_read_id;
					else
						mate_id = *b + max_read_id;
					if (data_tag[mate_id] > -2) {
						candi = data[*b].substr(i);
						candi_id = *b;
						break;
					}
					read_int = get_read_int_type(data[*b]);
					if (abs((int)read_hash[read_int].size() - ref_cov) < current_cov) {
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						candi = data[*b].substr(i);
						candi_id = *b;
					}
				}
				++b;				
			}//while(*b)
			if (candi_id > -1) {
				contig = contig + candi;
				read_int = get_read_int_type(data[candi_id]);
				candi_set = read_hash[read_int];
				for (int j = 0; j < candi_set.size(); ++j) {
					data_tag[candi_set[j]] = color_tag;
				}
				break;
			} 		
		}//for i
		if (candi_id == -1)
			break;
	} //while
*/
	return false;
}



string SplicingGraph::reverse_extend(ReadHash& read_hash, KmerHash& kmer_hash, const string& seed_contig, int color_tag) {

//cout << "reverse extend..." << endl;
	if (seed_contig.length() < g_read_length)
		return seed_contig;

	string contig = seed_contig;
//*
	int mate_id, r_d, candi_id,s_c;
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	vector<int> seed_readset;
	vector<int> tag_readset;
	kmer_int_type tag_kmer,kmer_int;

	while (1) {
		
		const string& reverse_kmer = contig.substr(0,g_kmer_length);
		kmer_int = kmer_to_int(reverse_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(0, g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		candi_id = -1;
		string::size_type start;
		s_c = g_read_length;
		if (g_is_paired_end) {
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				if (r_d >= max_read_id)
					mate_id = r_d - max_read_id;
				else
					mate_id = r_d + max_read_id;
				if (data_tag[mate_id] > -2) {
					start = data[r_d].find(reverse_kmer);
					if (start != string::npos && start < s_c && data[r_d].substr(start) == contig.substr(0, g_read_length-start)) {
						s_c = start;
						candi_id = r_d;
						if (s_c == 1)
							break;
						
					}
				}
			}
		}
		if (candi_id > -1) {
			contig = data[candi_id].substr(0, s_c) + contig;
			read_int = get_read_int_type(data[candi_id]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
			continue;
		}
		}
		s_c = g_read_length;
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				start = data[r_d].find(reverse_kmer);
				if (start != string::npos && data[r_d].substr(start) == contig.substr(0, g_read_length-start)) {
					read_int = get_read_int_type(data[r_d]);
					if (start < s_c) {
						s_c = start;
						candi_id = r_d;
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
					} else {
						if (start == s_c && abs((int)read_hash[read_int].size()-ref_cov) < current_cov) {
							s_c = start;
							candi_id = r_d;
							current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						}
					}	
						
				}
				
			}
		}
		if (candi_id > -1) {
			contig = data[candi_id].substr(0, s_c) + contig;
			read_int = get_read_int_type(data[candi_id]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
		} else {
			break;
		}


/*
	int mate_id, r_d, candi_id, candi_id1, s_c, s_c1;
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	vector<int> seed_readset;
	vector<int> tag_readset;
	kmer_int_type tag_kmer,kmer_int;

	while (1) {

		const string& reverse_kmer = contig.substr(0,g_kmer_length);
		kmer_int = kmer_to_int(reverse_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(0, g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		candi_id1 = -1;
		candi_id = -1;
		string::size_type start;
		s_c = g_read_length;
		s_c1 = g_read_length;
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				start = data[r_d].find(reverse_kmer);
				if (start == string::npos)
					continue;
				
				if (data[r_d].substr(start) == contig.substr(0, g_read_length-start)) {
					if (r_d < max_read_id && start < s_c1 && data_tag[r_d+max_read_id] > -2 ) {
						s_c1 = start;
						candi_id1 = r_d;
						if (s_c1 == 1)
							break;
					}
					read_int = get_read_int_type(data[r_d]);
					if (start < s_c) {
						s_c = start;
						candi_id = r_d;
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
					} else {
						if (start == s_c && abs((int)read_hash[read_int].size()-ref_cov) < current_cov) {
							s_c = start;
							candi_id = r_d;
							current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						}
					}	
						
				}
				
			}
		}
		if (candi_id1 > -1) {
			contig = data[candi_id1].substr(0, s_c1) + contig;
			read_int = get_read_int_type(data[candi_id1]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
		} else {
			if (candi_id > -1) {
				contig = data[candi_id].substr(0, s_c) + contig;
				read_int = get_read_int_type(data[candi_id]);
				candi_set = read_hash[read_int];
				for (int j = 0; j < candi_set.size(); ++j) 
				data_tag[candi_set[j]] = color_tag;
			} else {
				break;
			}
		}

/*
	string find_seq;
	int mate_id, r_d, candi_id;
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	vector<int> seed_readset;
	vector<int> tag_readset;
	kmer_int_type tag_kmer,kmer_int;

	while (1) {

		const string& reverse_kmer = contig.substr(0,g_kmer_length);
		kmer_int = kmer_to_int(reverse_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(0, g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		string candi;
		candi_id = -1;

		for (int i = 0; i < seed_readset.size(); ++i)
			data_find[seed_readset[i]] = true;

		for (int i = g_max_same_len-g_kmer_length; i >=g_min_same_len-g_kmer_length; --i) {
			const string& second_kmer = contig.substr(i,g_kmer_length);		
			tag_kmer = kmer_to_int(second_kmer);
	        	tag_readset = kmer_hash[tag_kmer];
			const string& ref_seq = contig.substr(0,i+g_kmer_length);
			vector<int>::iterator b = tag_readset.begin();
	        while (b != tag_readset.end()) {
		  		find_seq = data[*b].substr(g_read_length-i-g_kmer_length);
		  		
				if (data_find[*b] && data_tag[*b] < -2 && find_seq == ref_seq) {
					if (*b >= max_read_id)
						mate_id = *b - max_read_id;
					else
						mate_id = *b + max_read_id;
					if (data_tag[mate_id] > -2) {
						candi = data[*b].substr(0, g_read_length-i-g_kmer_length);
						candi_id = *b;
						break;
					}
					read_int = get_read_int_type(data[*b]);
					if (abs((int)read_hash[read_int].size() - ref_cov) < current_cov) {
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						candi = data[*b].substr(0, g_read_length-i-g_kmer_length);
						candi_id = *b;
					}
				} 		
				++b; 
			} //while

			if (candi_id > -1) {
				contig = candi + contig;
				read_int = get_read_int_type(data[candi_id]);
				candi_set = read_hash[read_int];
				for (int j = 0; j < candi_set.size(); ++j)
					data_tag[candi_set[j]] = color_tag;
				break;	
			} 		
		} //for (int i = g_max_same_len-g_kmer_length; i >=g_min_same_len-g_kmer_length; --i)
		
		if (candi_id == -1)
			break;
*/
	} //while (!stop_tag)
	
	return contig;
}


bool SplicingGraph::reverse_extend(ReadHash& read_hash, KmerHash& kmer_hash,  string& contig, const string& stop_seq, int color_tag) {

//cout << "reverse extend..." << endl;
	if (contig.length() < g_read_length)
		return false;
//*
	int mate_id, r_d, candi_id;
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	vector<int> seed_readset;
	vector<int> tag_readset;
	kmer_int_type tag_kmer,kmer_int;

	while (1) {
		const string& reverse_kmer = contig.substr(0,g_kmer_length);
		kmer_int = kmer_to_int(reverse_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(0, g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		candi_id = -1;
		string::size_type start;
		int s_c = g_read_length;
		if (g_is_paired_end) {
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				if (r_d >= max_read_id)
					mate_id = r_d - max_read_id;
				else
					mate_id = r_d + max_read_id;
				if (data_tag[mate_id] > -2) {
					start = data[r_d].find(reverse_kmer);
					if (start != string::npos && data[r_d].substr(start) == contig.substr(0, g_read_length-start)) {
						if (start < s_c ) {
							s_c = start;
							candi_id = r_d;
							if (s_c == 1)
								break;
						}
						int compatible_len = get_compatible_len(stop_seq, data[r_d].substr(0,g_min_same_len-1));
						if (compatible_len != 0){
							contig = data[r_d].substr(compatible_len) + contig.substr(g_read_length-start);
							data_tag[r_d] = color_tag;
							return true;
						}
					}
				}
			}
		}
		if (candi_id > -1) {
			contig = data[candi_id].substr(0, s_c) + contig;
			read_int = get_read_int_type(data[candi_id]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
			continue;
		}
		}
		s_c = g_read_length;
		for (int i = 0; i < seed_readset.size(); ++i) {
			r_d = seed_readset[i];
			if (data_tag[r_d] < -2) {
				start = data[r_d].find(reverse_kmer);
				if (start != string::npos && data[r_d].substr(start) == contig.substr(0, g_read_length-start)) {
					int compatible_len = get_compatible_len(stop_seq, data[r_d].substr(0,g_min_same_len-1));
					if (compatible_len != 0){
						contig = data[r_d].substr(compatible_len) + contig.substr(g_read_length-start);
						data_tag[r_d] = color_tag;
						return true;
					}
					read_int = get_read_int_type(data[r_d]);
					if (start < s_c) {
						s_c = start;
						candi_id = r_d;
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
					} else {
						if (start == s_c && abs((int)read_hash[read_int].size()-ref_cov) < current_cov) {
							s_c = start;
							candi_id = r_d;
							current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						}
					}	
						
				}
				
			}
		}
		if (candi_id > -1) {
			contig = data[candi_id].substr(0, s_c) + contig;
			read_int = get_read_int_type(data[candi_id]);
			candi_set = read_hash[read_int];
			for (int j = 0; j < candi_set.size(); ++j) {
				data_tag[candi_set[j]] = color_tag;
			}
		} else {
			break;
		}

/*
	string find_seq;
	int mate_id, r_d, candi_id;
	read_int_type read_int;
	int current_cov = 0;
	int ref_cov = 0;
	vector<int> candi_set;
	vector<int> seed_readset;
	vector<int> tag_readset;
	kmer_int_type tag_kmer,kmer_int;
	string::size_type start, start1;
	while (1) {

		const string& reverse_kmer = contig.substr(0,g_kmer_length);
		kmer_int = kmer_to_int(reverse_kmer);
		seed_readset = kmer_hash[kmer_int];
		const string& end_read = contig.substr(0, g_read_length);
		read_int = get_read_int_type(end_read);
		ref_cov = read_hash[read_int].size();
		current_cov = ref_cov*100000;
		string candi;
		candi_id = -1;

		for (int i = 0; i < seed_readset.size(); ++i)
			data_find[seed_readset[i]] = true;

		for (int i = g_max_same_len-g_kmer_length; i >=g_min_same_len-g_kmer_length; --i) {
			const string& second_kmer = contig.substr(i,g_kmer_length);		
			tag_kmer = kmer_to_int(second_kmer);
	        	tag_readset = kmer_hash[tag_kmer];
			const string& ref_seq = contig.substr(0,i+g_kmer_length);
			vector<int>::iterator b = tag_readset.begin();
	        while (b != tag_readset.end()) {
		  		find_seq = data[*b].substr(g_read_length-i-g_kmer_length);
		  		
				if (data_find[*b] && data_tag[*b] < -2 && find_seq == ref_seq) {

					int compatible_len = get_compatible_len(stop_seq, data[*b].substr(0,g_min_same_len-1));
					if (compatible_len != 0){
						contig = data[*b].substr(compatible_len) + contig.substr(i+g_kmer_length);
						data_tag[*b] = color_tag;
						return true;
					}
					if (*b >= max_read_id)
						mate_id = *b - max_read_id;
					else
						mate_id = *b + max_read_id;
					if (data_tag[mate_id] > -2) {
						candi = data[*b].substr(0, g_read_length-i-g_kmer_length);
						candi_id = *b;
						break;
					}
					read_int = get_read_int_type(data[*b]);
					if (abs((int)read_hash[read_int].size() - ref_cov) < current_cov) {
						current_cov = abs((int)read_hash[read_int].size() - ref_cov);
						candi = data[*b].substr(0, g_read_length-i-g_kmer_length);
						candi_id = *b;
					}
				} 		
				++b; 
			} //while

			if (candi_id > -1) {
				contig = candi + contig;
				read_int = get_read_int_type(data[candi_id]);
				candi_set = read_hash[read_int];
				for (int j = 0; j < candi_set.size(); ++j)
					data_tag[candi_set[j]] = color_tag;
				break;	
			} 		
		} //for (int i = g_max_same_len-g_kmer_length; i >=g_min_same_len-g_kmer_length; --i)
		
		if (candi_id == -1)
			break;
*/
	} //while (!stop_tag)
	
	return false;
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


void SplicingGraph::trim_graph(ReadHash& read_hash) {

	set_coverage_of_nodes(read_hash);
	bool has_trimed = false;
	map<pair_t, double> edge_coverages;
	vector<pair_t> edges;
	get_coverage_of_edges(read_hash, edge_coverages, edges);
	map<node_id, double> total_out_coverages;
	map<node_id, double> total_in_coverages;


	for (int i = 0; i < edges.size(); ++i) {

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

	for (int i = 0; i < edges.size(); ++i) {

		int source = edges[i].first;
		int target = edges[i].second;

		//if ((node_set[source].parents.empty() && node_set[source].children.size() == 1) || 
		//	(node_set[target].children.empty() && node_set[target].parents.size() == 1)) {

				const string& check_edge = get_edge_sequence(source, target);
				if ((int)check_edge.length() < g_read_length)
					continue;
				double e_cov = edge_coverages[edges[i]];
				double flanking_node_cov = node_set[source].node_coverage > node_set[target].node_coverage ? node_set[source].node_coverage : node_set[target].node_coverage;
				
				if (e_cov < g_min_junction_coverage ||
					e_cov < g_min_ratio_welds * flanking_node_cov ||
					e_cov < g_min_ratio_branch * total_out_coverages[source] ||
					e_cov < g_min_ratio_branch * total_in_coverages[target] ||
					(total_in_coverages.find(source) != total_in_coverages.end() && e_cov < g_min_ratio_in_out * total_in_coverages[source]) ||
					(total_out_coverages.find(target) != total_out_coverages.end() && e_cov < g_min_ratio_in_out * total_out_coverages[target])
					) {

						has_trimed = true;
						set_reads_tag(read_hash, check_edge, -4);
						node_set[source].delete_child(target);
						node_set[target].delete_parent(source);
				}

		//}

	}
	set_parents();
}



void SplicingGraph::set_coverage_of_nodes(ReadHash& read_hash) {

	for (int i = 0; i < node_sum; ++i) {

		if (node_set[i].sequence.length() < g_read_length) {
			node_set[i].node_coverage = 0;
			continue;
		}
		//set_reads_pos_in_node(read_hash, i);
		node_set[i].node_coverage = compute_coverage(read_hash, node_set[i].sequence);

	}

}


void SplicingGraph::get_coverage_of_edges(ReadHash& read_hash, map<pair_t, double>& edge_coverages, vector<pair_t>& edges) {

	edges.clear();
	edge_coverages.clear();

	for (int i = 0; i < node_sum; ++i) {

		if (node_set[i].children.empty())
			continue;

		for (int j = 0; j < (int)node_set[i].children.size(); ++j)
			edges.push_back(pair_t(i, node_set[i].children[j]));

	}

	for (int i = 0; i < (int)edges.size(); ++i) {

		int source = edges[i].first;
		int target = edges[i].second;
		const string& edge = get_edge_sequence(source, target);
		edge_coverages[edges[i]] = compute_coverage(read_hash,edge);

	}

}


void SplicingGraph::get_coverage_of_edges(ReadHash& read_hash, map<pair_t, double>& edge_coverages) {

	vector<pair_t> edges;

	for (int i = 0; i < node_sum; ++i) {

		if (node_set[i].children.empty())
			continue;

		for (int j = 0; j < (int)node_set[i].children.size(); ++j)
			edges.push_back(pair_t(i, node_set[i].children[j]));

	}

	for (int i = 0; i < (int)edges.size(); ++i) {

		int source = edges[i].first;
		int target = edges[i].second;
		const string& edge = get_edge_sequence(source, target);
		edge_coverages[edges[i]] = compute_coverage(read_hash, edge);

	}

}


float SplicingGraph::compute_coverage(ReadHash& read_hash, const string& sequence) {


	if (sequence.length() < g_read_length)
		return 0.0;
	float total = 0.0;
/*
	const int cov_size = static_cast<int>(sequence.length()) - g_read_length + 1;
	vector<int> seq_cov;
	seq_cov.reserve(cov_size);
	read_int_type read_int;
	for (int i = 0; i <= sequence.length()-g_read_length; ++i) {
		const string real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		seq_cov.push_back(read_hash[read_int].size());		
	}
	sort(seq_cov.begin(), seq_cov.end());
	int quantile = static_cast<int>(0.05*seq_cov.size() + 0.5);
	vector<int>::iterator start = seq_cov.begin();
	vector<int>::iterator end = seq_cov.end();
	if (quantile > 0) {
		start = seq_cov.begin() + quantile;
		end = seq_cov.end() - quantile;
	}
	for (; start != end; ++start)
		total = total + *start;

	total = total / (seq_cov.size() - 2*quantile);
	return total;
*/
//*
	vector<int> read_set;
	for (int i = 0; i <= sequence.length()-g_read_length; ++i) {
		const string& read = sequence.substr(i, g_read_length);
		read_int_type read_int =  get_read_int_type(read);
		if (read_hash[read_int].size() > 0) {
			read_set=read_hash[read_int];
			for (int j = 0; j < read_set.size(); ++j) {
				if (data_tag[read_set[j]] == -2)
					continue;
				total = total + 1.0;
			}
		}
	}

	total = total * g_read_length / sequence.length();
	return total;
	//*/
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

	return edge;

}


void SplicingGraph::get_transcripts(ReadHash& read_hash, vector<pair<string,float> >& transcripts) {
	cout << "get_transcripts..." << endl;
	set_coverage_of_nodes(read_hash);
	TreeStruct isotree;
	TreeNode root_node(node_order[node_order.size()-1],0);
 	int root_id = isotree.add_node(root_node);	
	vector<pair<node_id,float> > leaf;
	vector<int> tree_leaf_id;
	vector<vector<int> > trans_node_vec;
	vector<int> dealed;

	for (int i = 0; i < node_set[node_order[node_order.size()-1]].children.size(); ++i) {

		node_id root_child = node_set[node_order[node_order.size()-1]].children[i];

		if (node_set[root_child].node_coverage > g_min_ratio_non_error) {
		
			float node_out = 0;
			for (int j = 0; j < node_set[root_child].children.size(); ++j) {
				const string& edge_out = get_edge_sequence(root_child, node_set[root_child].children[j]);
				node_out = node_out + compute_coverage(read_hash, edge_out);
			}

			TreeNode leaf_node(root_child, node_out);	
			if (g_is_paired_end) {
			for (int j = 0; j <node_set[root_child].cov_reads.size(); ++j) {
				int mate_id = node_set[root_child].cov_reads[j].first;
				int rid = -1;
				if (mate_id >= max_read_id)
					rid = mate_id - max_read_id;
				else
					rid = mate_id + max_read_id;
				if ( data_tag[rid] >= 0 && data_tag[rid] != root_child ) {
					leaf_node.add_candi(data_tag[rid], 1);//////
				}
			}
			}
			int leaf_id = isotree.add_node(leaf_node);		
			isotree.add_child(root_id, leaf_id);
			pair <node_id, float> leaf_info;
			leaf_info.first = root_child;	
			leaf_info.second = node_out;	
			leaf.push_back(leaf_info);	
			tree_leaf_id.push_back(leaf_id);
		} else {
			for (int j = 0; j <node_set[root_child].cov_reads.size(); ++j) {
				int rid = node_set[root_child].cov_reads[j].first; 
				data_tag[rid] = -4;
			}
		}

	}

	for (int i = 1; i < (int)node_order.size()-1; ++i) {

		node_id x = node_order[node_order.size()-1-i];
/*
		if (node_set[x].node_coverage <= g_min_ratio_non_error ){
			for (int j = 0; j <node_set[x].cov_reads.size(); ++j)
				data_tag[node_set[x].cov_reads[j].first] = -1;
			continue;
		}
*/
		dealed.push_back(x);

		float total_out = 0.0;
		vector<pair<node_id, float> > out_edge;
		bool has_t = false;

		for (int j = 0; j < node_set[x].children.size(); ++j) {

			const string& edge = get_edge_sequence(x, node_set[x].children[j]);
			float edge_cov = compute_coverage(read_hash, edge);
			
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
		for (int j = 0; j < leaf.size(); ++j) {

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

		int out_sum = out_edge.size();
		int in_sum = in_edge.size();
		vector<float> adj(out_sum*in_sum, 0);

		solve_milp(out_edge, in_edge, adj, in_candi, in_priority);

		for (int j = 0; j < out_sum; ++j) {

			for (int k =0; k < in_sum; k++) {

				if (adj[k*out_sum+j] > 0) {				
					TreeNode c_node(out_edge[j].first, adj[k*out_sum+j]);
					c_node.candi = isotree.treenode_set[tree_leaf_id[in_edge[k].first]].candi;
					c_node.priority = isotree.treenode_set[tree_leaf_id[in_edge[k].first]].priority;
					c_node.dele_candi(out_edge[j].first);
					if (g_is_paired_end){
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

		
		for (int j = 0 ; j < in_sum; ++j) {

			for (int k = 0; k <in_sum-j-1; k++) {

				if (in_edge[k].first < in_edge[k+1].first) {

					pair<int, float> temp;
					temp = in_edge[k];
					in_edge[k] = in_edge[k+1];
					in_edge[k+1] = temp;
				
				}

			}

		}

		for(int j = 0; j < in_sum; ++j) {

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
                
		for (int i = 0; i < w_trans.size(); ++i) {
		
			string tran_seq;
			int exon_sum = 0;
			for (int j = 0; j < trans_vec[i].size(); ++j) {

				if (trans_vec[i][j] == node_order[0] || trans_vec[i][j] == node_order[node_order.size()-1])
					continue;


				if (j > 0) {				
					tran_seq = tran_seq + node_set[trans_vec[i][j]].sequence;
					++exon_sum;//cout << trans_vec[i][j] << ", "; ////
				}
				else {				
					tran_seq = node_set[trans_vec[i][j]].sequence;
					++exon_sum;//cout << trans_vec[i][j] << ", "; ////
				}

		
			}
			if (w_trans[i] < g_min_trans_cov || tran_seq.length() < 150+50*exon_sum) {
				for (int j = 0; j < trans_vec[i].size(); ++j) {
					for (int it = 0; it < node_set[trans_vec[i][j]].cov_reads.size(); it++)
						data_tag[node_set[trans_vec[i][j]].cov_reads[it].first] = -4;	
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
	check_transcripts(read_hash, transcripts);
}



void SplicingGraph::get_transcripts1(ReadHash& read_hash, vector<pair<string,float> >& transcripts) {
	cout << "get_transcripts..." << endl;
	set_coverage_of_nodes(read_hash);
	TreeStruct isotree;
	TreeNode root_node(node_order[node_order.size()-1],0);
 	int root_id = isotree.add_node(root_node);	
	vector<pair<node_id,float> > leaf;
	vector<int> tree_leaf_id;
	vector<vector<int> > trans_node_vec;
	vector<int> dealed;

	for (int i = 0; i < node_set[node_order[node_order.size()-1]].children.size(); ++i) {

		node_id root_child = node_set[node_order[node_order.size()-1]].children[i];

		if (node_set[root_child].node_coverage > g_min_ratio_non_error) {
		
			float node_in = 0;
			float node_out = 0;

			for (int j = 0; j < node_set[root_child].parents.size(); ++j) {
				const string& edge_in = get_edge_sequence(node_set[root_child].parents[j], root_child);
				node_in = node_in + compute_coverage(read_hash, edge_in);
			}

			for (int j = 0; j < node_set[root_child].children.size(); ++j) {
				const string& edge_out = get_edge_sequence(root_child, node_set[root_child].children[j]);
				node_out = node_out + compute_coverage(read_hash, edge_out);
			}


			TreeNode leaf_node(root_child, node_out-node_in);	

			
			for (int j = 0; j <node_set[root_child].cov_reads.size(); ++j) {
				int mate_id = node_set[root_child].cov_reads[j].first;
				int rid;
				if (mate_id >= max_read_id)
					rid = mate_id - max_read_id;
				else
					rid = mate_id + max_read_id;
				if ( data_tag[rid] >= 0 && data_tag[rid] != root_child && (!node_set[root_child].is_child(data_tag[rid])) ) {
					leaf_node.add_candi(data_tag[rid], 1);//////
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
			for (int j = 0; j <node_set[root_child].cov_reads.size(); ++j) {
				int rid = node_set[root_child].cov_reads[j].first; 
				data_tag[rid] = -4;
			}
		}

	}

	for (int i = 1; i < (int)node_order.size()-1; ++i) {

		node_id x = node_order[node_order.size()-1-i];

/*
		if (node_set[x].node_coverage <= g_min_ratio_non_error ){
			for (int j = 0; j <node_set[x].cov_reads.size(); ++j)
				data_tag[node_set[x].cov_reads[j].first] = -1;
			continue;
		}
*/
		dealed.push_back(x);

		float total_out = 0.0;
		vector<pair<node_id, float> > out_edge;
		bool has_t = false;

		for (int j = 0; j < node_set[x].children.size(); ++j) {

			const string& edge = get_edge_sequence(x, node_set[x].children[j]);
			float edge_cov = compute_coverage(read_hash, edge);
			
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
		for (int j = 0; j < leaf.size(); ++j) {

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

		for (int j = 0; j < out_sum; ++j) {

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

		
		for (int j = 0 ; j < in_sum; ++j) {

			for (int k = 0; k <in_sum-j-1; k++) {

				if (in_edge[k].first < in_edge[k+1].first) {

					pair<int, float> temp;
					temp = in_edge[k];
					in_edge[k] = in_edge[k+1];
					in_edge[k+1] = temp;
				
				}

			}

		}

		for(int j = 0; j < in_sum; ++j) {

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
                
		for (int i = 0; i < w_trans.size(); ++i) {
		
			string tran_seq;
		
			for (int j = 0; j < trans_vec[i].size(); ++j) {

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
				for (int j = 0; j < trans_vec[i].size(); ++j) {
					for (int it = 0; it < node_set[trans_vec[i][j]].cov_reads.size(); it++)
						data_tag[node_set[trans_vec[i][j]].cov_reads[it].first] = -4;	
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
	check_transcripts(read_hash, transcripts);
}



void SplicingGraph::topological_sort() {
	
	Node s,t;
	s.sequence = "start";
	t.sequence = "end";
	int s_id = add_node(s);
	int t_id = add_node(t);
	vector<int> node_color;

	for (int i = 0; i < node_sum; ++i) {

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

		for (int j = 0 ; j < node_set[i].children.size(); ++j) {

			if (node_color[node_set[i].children[j]] == 0)
				dfs_visit(node_set[i].children[j], node_color);

		}

		node_order.push_back(i);

	}

}



int SplicingGraph::get_total_amount_of_kmers() {

	int kmer_count = 0;

	for (unsigned int i = 0; i < node_sum; ++i)
		kmer_count += node_set[i].sequence.length();

	return (kmer_count - g_kmer_length);

}


/*
void SplicingGraph::solve_milp(vector<pair<node_id, float> >& out_edge, vector<pair<node_id, float> >& in_edge, vector<float>& adj, vector<vector<int> >& in_candi, vector<vector <float> >& in_priority) {
		
	int out_sum = out_edge.size();
	int in_sum = in_edge.size();
	if (out_sum == 1 && in_sum == 1){
		adj[0] = out_edge[0].second;
		return;
	}
	if (out_sum == 1) {
		for (int j = 0; j < in_sum; ++j)
			adj[j] = in_edge[j].second;
		return;
	}
	if (in_sum == 1) {
		for (int j = 0; j < out_sum; ++j)
			adj[j] = out_edge[j].second;
		return;
	}
		
	bool assinged = false;
	int in_t = 0;
	int out_t = 0;
	vector<bool> in_tag;
	vector<bool> out_tag;
	for (int i = 0; i < in_sum; ++i)
		in_tag.push_back(true);
	for (int i = 0; i < out_sum; ++i)
		out_tag.push_back(true);
		//新方法：
	bool break_while = true;
	while (break_while) {
		break_while = false;
		int y;
		for (int i = 0; i < out_sum; ++i) {
			float max_candi = 0;
			int candi_id = -1;
			for (int j = 0; j < in_sum; ++j) {
				for (int k = 0; k < in_candi[j].size(); k++) {
					if (in_candi[j][k] == out_edge[i].first) {
						if (in_priority[j][k] > max_candi && in_tag[j] == 0) {
							max_candi = in_priority[j][k];
							candi_id = j;
							y = k;
						}
						break;
					}
				}
			}
			if (max_candi > 0 && (in_tag[candi_id] || out_tag[i])) {
				if (out_edge[i].second > in_edge[candi_id].second) 					
					adj[out_sum*candi_id+i] = in_edge[candi_id].second;					
				 else
					adj[out_sum*candi_id+i] = out_edge[i].second;
				in_tag[candi_id] = false;
				out_tag[i] = false;
				in_priority[candi_id][y] = 0;
				break_while = true;
			} 
		}
	
	}	
cout << "1111" << endl;	
	int max_cov = -1;
	int max_id = -1;
	for (int i = 0; i < out_sum; ++i) {
		if (!out_tag[i])
			continue;
		for (int j = 0; j < in_sum; ++j) {
			if (in_tag[j]) {
				if (out_edge[i].second > in_edge[j].second){ 					
					adj[out_sum*j+i] = in_edge[j].second;	
					out_edge[i].second = out_edge[i].second - in_edge[j].second;
					in_edge[j].second = 0;
				}else {
					adj[out_sum*j+i] = out_edge[i].second;
					in_edge[j].second = in_edge[j].second - out_edge[i].second;
					out_edge[i].second = 0;
				}
				in_tag[j] = false;
				out_tag[i] = false;
				break;
			}
		}
		cout << "a" << endl;
		cout << in_sum << endl;
		cout << out_sum << endl;
		if (out_tag[i]) {
			max_cov = -1;
			for (int j = 0; j < in_sum; ++j) {
				if (in_edge[j].second > max_cov) {
					max_cov = in_edge[j].second;
					max_id = j;
				}
			}
			cout << max_cov <<" " << max_id << endl;
			if (out_edge[i].second > in_edge[max_id].second){ 					
				adj[out_sum*max_id+i] = in_edge[max_id].second;				
				out_edge[i].second = out_edge[i].second - in_edge[max_id].second;
				in_edge[max_id].second = 0;
			}else {
				adj[out_sum*max_id+i] = out_edge[i].second;
				in_edge[max_id].second = in_edge[max_id].second - out_edge[i].second;
				out_edge[i].second = 0;
			}
			in_tag[max_id] = false;
			out_tag[i] = false;
		}
		cout << "b" << endl;
	}

cout << "222" << endl;		
	for (int i = 0; i < in_sum; ++i) {
		if (!in_tag[i])
			continue;
		for (int j = 0; j < out_sum; ++j) {
			if (out_tag[j]) {
				if (in_edge[i].second > out_edge[j].second){ 					
					adj[out_sum*i+j] = out_edge[j].second;	
					in_edge[i].second = in_edge[i].second - out_edge[j].second;
					out_edge[j].second = 0;
				}else {
					adj[out_sum*i+j] = in_edge[i].second;
					out_edge[j].second = out_edge[j].second - in_edge[i].second;
					in_edge[i].second = 0;
				}
				in_tag[i] = false;
				out_tag[j] = false;
				break;
			}
		}
		if (in_tag[i]) {
			max_cov = -1;
			for (int j = 0; j < out_sum; ++j) {
				if (out_edge[j].second > max_cov) {
					max_cov = out_edge[j].second;
					max_id = j;
				}
			}
			if (in_edge[i].second > out_edge[max_id].second){ 					
				adj[out_sum*i+max_id] = out_edge[max_id].second;				
				in_edge[i].second = in_edge[i].second - out_edge[max_id].second;
				out_edge[max_id].second = 0;
			}else {
				adj[out_sum*i+max_id] = in_edge[i].second;
				out_edge[max_id].second = out_edge[max_id].second - in_edge[i].second;
				in_edge[i].second = 0;
			}
			out_tag[max_id] = false;
			in_tag[i] = false;
		}
	}
cout << "333" << endl;		
}
*/


void SplicingGraph::solve_milp(vector<pair<node_id, float> >& out_edge, vector<pair<node_id, float> >& in_edge, vector<float>& adj, vector<vector<int> >& in_candi, vector<vector <float> >& in_priority) {
		
	int out_sum = out_edge.size();
	int in_sum = in_edge.size();
	if (out_sum == 1 || in_sum == 1) {
		if (out_sum == 1 && in_sum == 1)			
			adj[0] = out_edge[0].second;
		else {
			if (out_sum == 1) {
				for (int j = 0; j < in_sum; ++j)
					adj[j] = in_edge[j].second;
			}
			if (in_sum == 1) {
				for (int j = 0; j < out_sum; ++j)
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
	while(1){
		int y = -1;
		bool stop_tag = true;
		for (int i = 0; i < out_sum; ++i) {
			float max_candi = 0;
			int candi_id = -1;
			for (int j = 0; j < in_sum; ++j) {

				for (int k = 0; k < in_candi[j].size(); k++) {
					if (in_candi[j][k] == out_edge[i].first) {
						if (in_priority[j][k] > max_candi && in_tag[j] == 0) {
							max_candi = in_priority[j][k];
							candi_id = j;
							y = k;
						}
						break;
					}
				}
			}

			//if (max_candi > 1) {

			if (max_candi > 1 && (in_tag[candi_id] == 0 || out_tag[i] == 0)) {

				if (out_edge[i].second > in_edge[candi_id].second)					
					adj[out_sum*candi_id+i] = in_edge[candi_id].second;
				else
					adj[out_sum*candi_id+i] = out_edge[i].second;
				in_priority[candi_id][y] = 0;
				in_tag[candi_id] =1;
				out_tag[i] = 1;
				stop_tag = false;
			}

		}
		if(stop_tag)
			break;
	}
		for (int i = 0; i < out_sum; ++i) {

			if (out_tag[i] == 0) {
				bool assign_tag = false;
				for (int j = 0; j < in_sum; ++j) {
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
					for (int j = 0; j < in_sum; ++j) {
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


		
		for (int i = 0; i < in_sum; ++i) {

			if (in_tag[i] == 0) {
				bool assign_tag = false;
				for (int j = 0; j < out_sum; ++j) {
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
					for (int j = 0; j < out_sum; ++j) {
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
		for (int i=0; i< out_sum; ++i)
			sum_vec.push_back(out_edge[i].second);

		for (int i = 0; i < out_sum; ++i) {
			float re_sum = 0.0;
			for (int j = 0; j < in_sum; ++j) 
				re_sum = re_sum + adj[j*out_sum+i];
			for (int j = 0; j < in_sum; ++j) {
				adj[j*out_sum+i] = sum_vec[i]*adj[j*out_sum+i] / re_sum;
			}
		}
	}
		
}




void SplicingGraph::recover_reads() {

	for (int i = 0; i < node_sum; ++i) {
		
		for (int it = 0; it < node_set[i].cov_reads.size(); it++) {
			data_tag[node_set[i].cov_reads[it].first] = -4;
		}
	}



}

int SplicingGraph::get_read_id(node_id p, int pos){


	for (int i = pos; i < node_set[p].sequence.length()-g_read_length; ++i) {

		const string& read = node_set[p].sequence.substr(i, g_read_length);
		for(int it = 0; it < node_set[p].cov_reads.size(); it++){
			int rid = node_set[p].cov_reads[it].first;
			if (data[rid] == read)
				return rid;
		}
	}

	return -1;
}


void SplicingGraph::get_edge_reads(ReadHash & read_hash, node_id source, node_id target, vector<int>& edge_reads) {

	if (source < 0 || source >= node_sum)
		return;
	if (target < 0 || target >= node_sum)
		return;
	const string& sequence = get_edge_sequence(source, target);
	read_int_type read_int;
	vector<int> read_set;
	for (int i = 0; i <= sequence.length()-g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_set = read_hash[read_int];
		edge_reads.insert(edge_reads.end(), read_set.begin(), read_set.end());
	}
}



void SplicingGraph::set_reads_tag(ReadHash& read_hash, const string& sequence, int tag) {

	if (sequence.length() < g_read_length)
		return;
	read_int_type read_int;
	vector<int> read_set;
	for (int i = 0; i <= sequence.length() - g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_set = read_hash[read_int];
		for (int j = 0; j < read_set.size(); ++j) {
			if (data_tag[read_set[j]] == -2)
				continue;
			data_tag[read_set[j]] = tag;
		}
	}
}


void SplicingGraph::set_reads_tag(ReadHash& read_hash, const string& sequence, vector<int>& map_reads, int tag) {

	if (sequence.length() < g_read_length)
		return;

	read_int_type read_int;
	vector<int> read_set;
	for (int i = 0; i <= sequence.length() - g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_set = read_hash[read_int];
		map_reads.insert(map_reads.end(), read_set.begin(), read_set.end());
		for (int j = 0; j < read_set.size(); ++j) {
			if (data_tag[read_set[j]] == -2)
				continue;
			data_tag[read_set[j]] = tag;
		}
	}
}



void SplicingGraph::map_reads_to_sequence(ReadHash& read_hash, const string& sequence, vector<pair<int, int> >& reads_pos) {

	if (sequence.length() < g_read_length)
		return;

	read_int_type read_int;
	vector<int> read_set;
	//read_set.reserve(1000);
	pair<int, int> r_p;
	for (int i = 0; i <= sequence.length() - g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_set = read_hash[read_int];
		r_p.second = i;
		for (int j = 0; j < read_set.size(); ++j) {
			if (data_tag[read_set[j]] ==-2)
				continue;
			r_p.first = read_set[j];
			data_tag[r_p.first] = i;
			reads_pos.push_back(r_p);
		}		
	}	
}



void SplicingGraph::check_transcripts(ReadHash& read_hash, vector<pair<string,float> >& transcripts) {

	cout << "check_transcripts.." << endl;

	vector<int> add_trans;
	vector<int> set_reads_vec;
	vector<int> recover_reads_vec;
	set<int> used_reads;
	float trans_cov;
///*
	string max_transcript = "N";
	int max_tran_id = -1;
	for (int i = 0; i < transcripts.size(); ++i) 
		if (transcripts[i].first.length() > max_transcript.length()) {
			max_transcript = transcripts[i].first;
			max_tran_id = i;
		}

//*/
	vector<pair<int, int> > reads_pos;
	//reads_pos.reserve(10000);
	for (int i = 0; i < transcripts.size(); ++i) {
		string transcript = transcripts[i].first;
		reads_pos.clear();
		map_reads_to_sequence(read_hash, transcript, reads_pos);
		bool is_add = true;
///*
		if (is_seq_similar(max_transcript, transcript, 'F',0.05) && i != max_tran_id)
			is_add = false;
		if (is_seq_similar(max_transcript, transcript,'R', 0.05) && i != max_tran_id)
			is_add = false;

		trans_cov = reads_pos.size()*g_read_length* 1.0 / transcript.length();
		if ((transcript.length() < 500 && trans_cov < g_min_trans_cov*4) || (transcript.length() < 800 && trans_cov < g_min_trans_cov*2) || (transcript.length() < 1000 && trans_cov < g_min_trans_cov*1.5))
			is_add = false;
		if ((!is_add) || trans_cov < g_min_trans_cov || transcript.length() < g_min_transcript_length) {

			for (int j = 0; j < reads_pos.size(); ++j)
				data_tag[reads_pos[j].first] = -4;
			continue;
		}

		if (g_is_paired_end) {
			
			float pair_end_sum = 0.0;		
			vector<pair<int, int> > maybe_used;
			for (int j = 0; j < reads_pos.size(); ++j) {

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
					data_tag[read_id] = -4;
				}

			}

			trans_cov = pair_end_sum * g_read_length / transcript.length();
cout << pair_end_sum  << " " << reads_pos.size() << endl;
			if ((transcript.length() < 500 && trans_cov < g_min_trans_cov*4) ||
				(transcript.length() < 700 && trans_cov < g_min_trans_cov*2) || 				(transcript.length() < 1000 && trans_cov < g_min_trans_cov*1.2))
				is_add = false;
			if (is_add && pair_end_sum > 10) {
				int total_cov_len = 0;
				int current_cov_len = g_read_length;
				int current_cov_start = maybe_used[0].second;
				bool add_tag = true;
				
				for (int j = 0; j < maybe_used.size()-1; ++j) {
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
cout << "total_cov_len " << total_cov_len << endl;
cout << transcript.length() << endl;
				if ((total_cov_len > 0.8*transcript.length() && transcript.length() >1200) || (total_cov_len>0.9*transcript.length() && transcript.length() > 500) || total_cov_len > 0.99*transcript.length() || total_cov_len > 1500){
					transcripts[i].second = pair_end_sum * g_read_length / transcript.length();
					//for (int j = 0; j < maybe_used.size(); ++j)
						//data_used[maybe_used[j].first] = true;
				} else {
					is_add = false;
				}
			} else { 
				is_add = false;
			}					
		}
		
		if (is_add) {
			add_trans.push_back(i);			
		}
		for (int j = 0; j < reads_pos.size(); ++j)
			data_tag[reads_pos[j].first] = -4;

	}

	vector<pair<string,float> > temp_trans;
	for (int i = 0; i < add_trans.size(); ++i) {
		int tid = add_trans[i];
		temp_trans.push_back(transcripts[tid]);
	}
	transcripts.clear();	
	transcripts = temp_trans;

}





void SplicingGraph::set_reads_pos_in_node(ReadHash& read_hash, node_id p) {

	if (p < 0 || p >= node_sum )
		return;
	if (node_set[p].sequence.length() < g_read_length)
		return;

	const string& sequence = node_set[p].sequence;
/*
	if (node_set[p].cov_reads.size() > 0) {
		for (int i = 0; i < node_set[p].cov_reads.size(); ++i)
			data_tag[node_set[p].cov_reads[i].first] = -4;
		node_set[p].cov_reads.clear();
	}
*/
	node_set[p].cov_reads.clear();
	read_int_type read_int;
	vector<int> read_set;
	//read_set.reserve(1000);
	pair<int, int> r_p;
	for (int i = 0; i <= sequence.length() - g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_set = read_hash[read_int];
		r_p.second = i;
		for (int j = 0; j < read_set.size(); ++j) {
			r_p.first = read_set[j];
			data_tag[r_p.first] = i;
			node_set[p].cov_reads.push_back(r_p);
		}		
	}

}



void SplicingGraph::branch_extend_by_coverage(ReadHash& read_hash, KmerHash& kmer_hash) {

	int current_node_sum = node_sum;
	for (int i = 0; i < current_node_sum; ++i)
		forward_extend_by_coverage(read_hash, kmer_hash, i);
	current_node_sum = node_sum;
	for (int i = 0; i < current_node_sum; ++i)
		reverse_extend_by_coverage(read_hash, kmer_hash, i);

}


void SplicingGraph::forward_extend_by_coverage(ReadHash& read_hash, KmerHash& kmer_hash, node_id p) {

	if (p < 0 || p >= node_sum)
		return;

	set_forward_check_pos(read_hash, p);
	
	//cout << "check forward branch..." << endl;

	pair<int, int> right_id;
	read_int_type read_int;
//cout << node_set[p].forward_check_pos.size() << endl;
	while (node_set[p].forward_check_pos.size() > 0) {

		pair<int, int> check_pos = node_set[p].forward_check_pos.front();
		node_set[p].forward_check_pos.pop_front();
		if (node_set[p].sequence.length() < check_pos.second + g_read_length)
			break;
		if (check_pos.first <= g_read_length)
			break;
//cout << check_pos.first << " " << check_pos.second << " " << node_set[p].sequence.length() << endl;

		for (int i = check_pos.first; i < check_pos.second; ++i) {
			const string& seed = node_set[p].sequence.substr(i-g_read_length, g_read_length);
			read_int = get_read_int_type(seed);
			if (read_hash[read_int].size() == 0)
				continue;
			right_id.first = -1;
			right_id.second = -1;
			const string& forward_seq = forward_extend(read_hash, kmer_hash, seed);
			set_reads_tag(read_hash, forward_seq, -4);

			if ((g_is_paired_end &&(!pair_support(read_hash, forward_seq,p))) || forward_seq.length() <= g_read_length + 10) 
				continue;
			const string& reach_seq = forward_seq.substr(forward_seq.length()-g_min_same_len);
			string::size_type start2 = node_set[p].sequence.substr(i).find(reach_seq);
			int right_node = -1;
			if (start2 == string::npos) {
				for (int k = 0; k < node_sum; ++k) {
					if (k ==p)
						continue;
					start2 = node_set[k].sequence.find(reach_seq);
					if (start2 != string::npos) {
						if (g_is_paired_end && (!pair_support(read_hash, forward_seq,k)))
							continue;
						right_id.first = k;
						right_id.second = start2;
						right_node = k;
						break;
					}
				}

			} else {
				right_id.first = p;
				right_node = p;
				right_id.second = start2;
			}


			if (right_id.second == -1 && forward_seq.length() < g_min_exon_length + g_read_length)
				continue;

			if (right_id.second == -1) {//add branch
				Node node2 = node_set[p].sequence.substr(0, i + g_read_length);
				int n2 = add_node(node2);
				node_set[p].sequence = node_set[p].sequence.substr(i+g_read_length);
				for (int j = 0; j < node_set[p].parents.size(); ++j) {
					int p_i = node_set[p].parents[j];
					node_set[p_i].delete_child(p);
					node_set[p_i].add_child(n2);
				}
				node_set[n2].parents = node_set[p].parents;
				node_set[p].parents.clear();
				node_set[p].add_parent(n2);
				node_set[n2].add_child(p);
				if (node_set[p].sequence.length() < g_min_exon_length && node_set[p].children.size() == 0) {
					set_reads_tag(read_hash, node_set[p].sequence, -4);
					node_set[p].sequence = forward_seq.substr(g_read_length);
					set_reads_tag(read_hash, forward_seq, p);
					forward_extend_by_coverage(read_hash, kmer_hash, p);
					return;
				} else {
					Node node1;
					node1.sequence = forward_seq.substr(g_read_length);
					int n1 = add_node(node1);
					node_set[n2].add_child(n1);
					node_set[n1].add_parent(n2);
					set_reads_tag(read_hash, forward_seq, n1);
					forward_extend_by_coverage(read_hash, kmer_hash, p);
					return;
				}
	
			} else { //add bubble
				bool fail_add = false;
				int length = forward_seq.length() - g_read_length - right_id.second;
				int start1 = i;
				string bubble_seq;
				if (length > 0)
					bubble_seq = forward_seq.substr(g_read_length, length);
				if ( p == right_node) {
					if ((!fail_add) && start2 <= start1)
						fail_add = true;
					int distance = 0;
					if (!fail_add)
						distance = start2 - start1 - g_read_length;
					if (distance <= 0)
						fail_add = true;
					if (distance + length <= 4)
						fail_add = true;
					if (!fail_add) {
						Node node1, node2, node3;
						node_id n1 = -1;
						node_id n2 = -1;
						node_id n3 = -1;
						if (length > 0) {
							node1.sequence = bubble_seq;
							n1 = add_node(node1);
							set_reads_tag(read_hash, forward_seq, n1);
						}
						node2.sequence = node_set[p].sequence.substr(0, start1+g_read_length);
						node3.sequence = node_set[p].sequence.substr(start1+g_read_length, distance);
						node_set[p].sequence = node_set[p].sequence.substr(start2);
						n2 = add_node(node2);
						n3 = add_node(node3);
						for (int j = 0; j < node_set[p].parents.size(); ++j) {
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
						forward_extend_by_coverage(read_hash, kmer_hash, n3);
						forward_extend_by_coverage(read_hash, kmer_hash, p);
						return;
					}

				} else { //p != right_node
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
							set_reads_tag(read_hash, forward_seq, n1);
						}
						if (start1+g_read_length < node_set[p].sequence.length()) {
							node2.sequence = node_set[p].sequence.substr(0, start1+g_read_length);
							n2 = add_node(node2);
							node_set[p].sequence = node_set[p].sequence.substr(start1+g_read_length);
							for (int j = 0; j < node_set[p].parents.size(); ++j) {
								int p_i = node_set[p].parents[j];
								node_set[p_i].delete_child(p);
								node_set[p_i].add_child(n2);
							}
							node_set[n2].parents = node_set[p].parents;
							node_set[p].parents.clear();
							node_set[n2].add_child(p);
							node_set[p].add_parent(n2);
							if (length > 0) {
								node_set[n2].add_child(n1);
								node_set[n1].add_parent(n2);
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
							node_set[right_node].sequence = node_set[right_node].sequence.substr(start2);
							for (int j = 0; j < node_set[right_node].parents.size(); ++j) {
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
						forward_extend_by_coverage(read_hash, kmer_hash, p);
						return;
					}
				}


				if (fail_add && forward_seq.length() > g_min_exon_length + g_read_length){
						
					Node node2 = node_set[p].sequence.substr(0, i + g_read_length);
					int n2 = add_node(node2);
					node_set[p].sequence = node_set[p].sequence.substr(i + g_read_length);
					for (int j = 0; j < node_set[p].parents.size(); ++j) {
						int p_i = node_set[p].parents[j];
						node_set[p_i].delete_child(p);
						node_set[p_i].add_child(n2);
					}
					node_set[n2].parents = node_set[p].parents;
					node_set[p].parents.clear();
					node_set[p].add_parent(n2);
					node_set[n2].add_child(p);
					if (node_set[p].sequence.length() < g_min_exon_length && node_set[p].children.size() == 0) {
						set_reads_tag(read_hash, node_set[p].sequence, -4);
						node_set[p].sequence = forward_seq.substr(g_read_length);
						forward_extend_by_coverage(read_hash, kmer_hash, p);
						return;
					} else {
						Node node1;
						node1.sequence = forward_seq.substr(g_read_length);
						int n1 = add_node(node1);
						set_reads_tag(read_hash, forward_seq, n1);
						node_set[n2].add_child(n1);
						node_set[n1].add_parent(n2);							
						forward_extend_by_coverage(read_hash, kmer_hash, p);
						return;
					}

				} //if (fail_add)
			}//add bubble

		} // for

	}

}



void SplicingGraph::reverse_extend_by_coverage(ReadHash& read_hash, KmerHash& kmer_hash, node_id p) {

	if (p < 0 || p >= node_sum)
		return;

	set_reverse_check_pos(read_hash, p);
//cout << "reverse  check " << p << endl;

	pair<int, int> right_id;
	read_int_type read_int;
	while (node_set[p].reverse_check_pos.size() > 0) {

		pair<int, int> check_pos = node_set[p].reverse_check_pos.front();
		node_set[p].reverse_check_pos.pop_front();
		if (node_set[p].sequence.length() <= check_pos.second + g_read_length)
			break;
		if (check_pos.first <= 0)
			break;
//cout << check_pos.first << " " << check_pos.second <<" "<< node_set[p].sequence.length() << endl;
		for (int i = check_pos.first; i < check_pos.second; ++i) {
			const string& seed = node_set[p].sequence.substr(i, g_read_length);
			read_int = get_read_int_type(seed);
			if (read_hash[read_int].size() == 0)
				continue;
			right_id.first = -1;
			right_id.second = -1;
			const string& reverse_seq = reverse_extend(read_hash, kmer_hash, seed);
			set_reads_tag(read_hash, reverse_seq, -4);

			string::size_type start2;
			int right_node = -1;
			const string& reach_seq = reverse_seq.substr(0, g_min_same_len);
			for (int k = 0; k < node_sum; ++k) {
				if (k ==p)
					continue;
				start2 = node_set[k].sequence.find(reach_seq);
				if (start2 != string::npos) {
					if (g_is_paired_end && (!pair_support(read_hash, reverse_seq,k)))
						continue;
					right_id.first = k;
					right_id.second = start2;
					right_node = k;
					break;
				}
			}

			if ((g_is_paired_end && (!pair_support(read_hash, reverse_seq, p))) || reverse_seq.length() <= g_read_length + 10) 
				continue;
			
			if (right_id.second == -1 && reverse_seq.length() < g_min_exon_length + g_read_length)
				continue;

			if (right_id.second == -1) {//add branch
				Node node2 = node_set[p].sequence.substr(0, i);
				int n2 = add_node(node2);
				node_set[p].sequence = node_set[p].sequence.substr(i);
				for (int j = 0; j < node_set[p].parents.size(); ++j) {
					int p_i = node_set[p].parents[j];
					node_set[p_i].delete_child(p);
					node_set[p_i].add_child(n2);
				}
				node_set[n2].parents = node_set[p].parents;
				node_set[p].parents.clear();
				node_set[p].add_parent(n2);
				node_set[n2].add_child(p);
				if (node_set[n2].sequence.length() < g_min_exon_length && node_set[n2].parents.size() == 0) {
					set_reads_tag(read_hash, node_set[n2].sequence, -4);
					node_set[n2].sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
					set_reads_tag(read_hash, reverse_seq, n2);
					reverse_extend_by_coverage(read_hash, kmer_hash, p);
					return;
				} else {
					Node node1;
					node1.sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
					int n1 = add_node(node1);
					node_set[n1].add_child(p);
					node_set[p].add_parent(n1);
					set_reads_tag(read_hash, reverse_seq, n1);
					reverse_extend_by_coverage(read_hash, kmer_hash, p);
					return;
				}
	
			} else { //add bubble
				bool fail_add = false;
				int length = reverse_seq.length() - g_read_length - right_id.second;
				int start1 = i;
				string bubble_seq;
				if (length > 0)
					bubble_seq = reverse_seq.substr(right_id.second, length);
				set<node_id> checked;
				if (right_node == 0 || is_circle(right_node, p, checked))
					fail_add = true;							
				if (!fail_add) {
					Node node1, node2, node3;
					node_id n1 = -1;
					node_id n2 = -1;
					node_id n3 = -1;
					if (length > 0) {
						node1.sequence = bubble_seq;
						n1 = add_node(node1);
						set_reads_tag(read_hash, reverse_seq, n1);
					}
					if (start1 > 0) {
						node2.sequence = node_set[p].sequence.substr(0, start1);
						n2 = add_node(node2);
						node_set[p].sequence = node_set[p].sequence.substr(start1);
						for (int j = 0; j < node_set[p].parents.size(); ++j) {
							int p_i = node_set[p].parents[j];
							node_set[p_i].delete_child(p);
							node_set[p_i].add_child(n2);
						}
						node_set[n2].parents = node_set[p].parents;
						node_set[p].parents.clear();
						node_set[n2].add_child(p);
						node_set[p].add_parent(n2);
						
					} 
					if (length > 0) {
						node_set[n1].add_child(p);
						node_set[p].add_parent(n1);
					}
					
					if (start2 + right_id.second < node_set[right_node].sequence.length()) {
						node3.sequence = node_set[right_node].sequence.substr(0, start2+right_id.second);
						n3 = add_node(node3);
						node_set[right_node].sequence = node_set[right_node].sequence.substr(start2 + right_id.second);
						for (int j = 0; j < node_set[right_node].parents.size(); ++j) {
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
						if (start2 + right_id.second < node_set[right_node].sequence.length()) {
							node_set[n1].add_parent(n3);
							node_set[n3].add_child(n1);
						} else {
							node_set[n1].add_parent(right_node);
							node_set[right_node].add_child(n1);
						}
					} else {
						if (start2 + right_id.second < node_set[right_node].sequence.length()) {
							node_set[n3].add_child(p);
							node_set[p].add_parent(n3);
						} else {
							node_set[right_node].add_child(p);
							node_set[p].add_parent(right_node);
						}
					}
					reverse_extend_by_coverage(read_hash, kmer_hash, p);
					return;
				}

				if (fail_add && reverse_seq.length() > g_min_exon_length + g_read_length){
					Node node1;
					int n1 = -1;
					if (i == 0) {
						if (node_set[p].parents.size() == 0) {
							node_set[p].sequence = reverse_seq + node_set[p].sequence.substr(g_read_length);
							set_reads_tag(read_hash, reverse_seq, p);
						} else {	
							node1.sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
							n1 = add_node(node1);
							set_reads_tag(read_hash, reverse_seq, n1);
							node_set[n1].add_child(p);
							node_set[p].add_parent(n1);							
						}
					} else {
						Node node2 = node_set[p].sequence.substr(0, i);
						int n2 = add_node(node2);
						node_set[p].sequence = node_set[p].sequence.substr(i);
						for (int j = 0; j < node_set[p].parents.size(); ++j) {
							int p_i = node_set[p].parents[j];
							node_set[p_i].delete_child(p);
							node_set[p_i].add_child(n2);
						}
						node_set[n2].parents = node_set[p].parents;
						node_set[p].parents.clear();
						node_set[p].add_parent(n2);
						node_set[n2].add_child(p);
						node1.sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
						n1 = add_node(node1);
						set_reads_tag(read_hash, reverse_seq, n1);
						node_set[n1].add_child(p);
						node_set[p].add_parent(n1);
					}							
					reverse_extend_by_coverage(read_hash, kmer_hash, p);
					return;

				} //if (fail_add)
			}//add bubble

		} // for

	}

}


void SplicingGraph::set_forward_check_pos(ReadHash& read_hash, node_id p) {

//cout << "set check" << endl;
	node_set[p].forward_check_pos.clear();
	if (node_set[p].sequence.length() <= g_read_length + 1)
		return;

	set_reads_pos_in_node(read_hash, p);

	if (node_set[p].cov_reads.size() < 10)
		return;

	int node_len = node_set[p].sequence.length();
	vector<float> node_cov(node_len, 0.0);
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i)
		for (int j = node_set[p].cov_reads[i].second; j < node_set[p].cov_reads[i].second+g_read_length; ++j)
			node_cov[j] = node_cov[j] + 1.0;
///*
	int start, len, end;
	pair<int, int> fp;
	len = 0;
	start = 0;
	end = 0;

	for (int i = g_read_length; i < (int)node_cov.size()-1; ++i) {

		if (node_cov[i] <= node_cov[i+1]){

			if (len >= 10 && node_cov[end] < 0.5*node_cov[start]) {

				fp.first = start-g_read_length;
				fp.second = start-g_read_length+10;
				node_set[p].forward_check_pos.push_back(fp);
			}
			start = 0;
			end = 0;
			len = 0;
		
		} else {

			if (start > 0) {
				++len;
				end = i+1;
			} else {
				start = i+1;
			}
		}
	}

	int id = -1;
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i) {
		id = node_set[p].cov_reads[i].first;
		data_tag[id] = p;
	}
//*/
/*
	bool first_forward = true;
	int forward_start = 0;
	int forward_end = 0;
	bool forward_tag = false;
	vector<pair<int, int> > candi_forward;

	for (int i = 1; i < (int)node_cov.size()-g_read_length - 1; ++i) {

		if (node_cov[i] < node_cov[i+g_min_same_len-1] * 0.4){

			if (forward_tag) {
				pair<int, int> fp;
				fp.first = forward_start;
				fp.second = forward_end;
				node_set[p].forward_check_pos.push_back(fp);
				first_forward = true;
				forward_tag = false;
			}
		
		} else {

			if (node_cov[i+g_min_same_len-1] < node_cov[i] * 0.4){

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
*/

}


void SplicingGraph::set_reverse_check_pos(ReadHash& read_hash, node_id p) {

//cout << "set check" << endl;
	node_set[p].reverse_check_pos.clear();

	if (node_set[p].sequence.length() <= g_read_length + 1)
		return;

	set_reads_pos_in_node(read_hash, p);

	if (node_set[p].cov_reads.size() < 10)
		return;

	int node_len = node_set[p].sequence.length();
	vector<float> node_cov(node_len, 0.0);
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i)
		for (int j = node_set[p].cov_reads[i].second; j < node_set[p].cov_reads[i].second+g_read_length; ++j)
			node_cov[j] = node_cov[j] + 1.0;
//*	
	int start, len, end;
	pair<int, int> fp;
	start = 0;
	end = 0;
	len = 0;
	for (int i = g_read_length; i < (int)node_cov.size()-1; ++i) {

		if (node_cov[i] >= node_cov[i+1]){
			if (len >= 10 && node_cov[end] > 0.5*node_cov[start]) {
				fp.first = end-10;
				fp.second = end;
				node_set[p].reverse_check_pos.push_back(fp);
			}
			start = 0;
			end = 0;
			len = 0;
		
		} else {

			if (start > 0) {
				++len;
				end = i+1;
			} else {
				start = i+1;
			}
		}
	}

	int id = -1;
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i) {
		id = node_set[p].cov_reads[i].first;
		data_tag[id] = p;
	}
//*/
/*

	bool first_reverse = true;
	int reverse_start = 0;
	int reverse_end = 0;
	bool reverse_tag = false;
	vector<pair<int, int> > candi_reverse;
	for (int i = 1; i < (int)node_cov.size()-g_read_length - 1; ++i) {
		if (node_cov[i] < node_cov[i+g_min_same_len-1] * 0.4){ 

			if (first_reverse) {
				reverse_start = i;
				first_reverse = false;
			}
			reverse_end = i + g_min_same_len-1;
			reverse_tag = true;
		
		} else {

			if (node_cov[i+g_min_same_len-1] < node_cov[i] * 0.4){// || node_cov[i] - node_cov[i+g_min_same_len-1] > 200) {

				if (reverse_tag) {
					pair<int,int> rp;
					rp.first = reverse_start;
					rp.second = reverse_end;
					node_set[p].reverse_check_pos.push_back(rp);
					first_reverse = true;
					reverse_tag = false;
				}

			}else {

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
*/
}


