
#ifndef SPLICINGGRAPH_H
#define SPLICINGGRAPH_H

//extern "C" {
//#include "LingoProgram.h"
//}

#include"KmerUtility.h"
#include"KmerHash.h"
#include "ReadHash.h"
#include "TreeStruct.h"
#include<algorithm>
#include<assert.h>
#include<fstream>
#include<map>
#include<numeric>
#include<iomanip>
#include<set>
#include<string>
#include <stdlib.h>
#include<sstream>
#include<vector>
#include<list>




using namespace std;


//typedef pair<int,int> pair_t;
//extern set<pair_t> Inhibit_edges;
//extern vector<vector<int>> Pair_edges;

typedef int node_id;


int get_compatible_len(const string& str1, const string& str2);
bool is_seq_similar(const string& str1, const string& str2, char mode = 'F', double sim_error = 0.35);
bool is_aligned(const string& str1, const string& str2, char mode = 'F', int tag = 0);
bool compatible(const string& str1, const string& str2);




class Node {

		//friend class  boost::serialization::access;
		//template<class Archive>
		//void serialize(Archive & ar, const unsigned int version) {
			//ar & sequence;
			//ar & parents; 
			//ar & children;
			//ar & node_coverage;
		//}

	public:

		string sequence;
		vector<node_id> parents;
		vector<node_id> children;
		float node_coverage;
		vector<pair<int, int> > cov_reads;
		list<pair<int, int> > forward_check_pos;
		list<pair<int, int> > reverse_check_pos;

		Node(): sequence("") {};

		Node(const string& mysequence){

			sequence = mysequence;

		}

		Node(const Node& node){

			sequence = node.sequence;
			parents = node.parents;
			children = node.children;
			node_coverage = node.node_coverage;
			cov_reads = node.cov_reads;
			forward_check_pos = node.forward_check_pos;
			reverse_check_pos = node.reverse_check_pos;
		}


		void set_sequence(const string& mysequence){

			sequence = mysequence;

		}

		void set_cov(float mynode_cov) {

			node_coverage = mynode_cov;

		}

		float get_node_coverage() {

			return node_coverage;

		}


		string get_sequence() {

			return sequence;

		}

		bool add_child(node_id child) {

			if(child < 0)
				return false;

			if(!children.empty()) {

				for (int i = 0; i < children.size(); i++) 
					if (children[i] == child)
						return false;
			}

			this -> children.push_back(child);
			return true;

		}

		bool add_parent(node_id parent) {

			if(parent < 0)
				return false;

			if(!parents.empty()) {

				for (int i = 0; i < parents.size(); i++)
					if (parents[i] == parent)
						return false;
			}

			this -> parents.push_back(parent);
			return true;

		}

		bool is_child(node_id child) {

			if (child < 0)
				return false;

			for(int i = 0; i < children.size(); i++) {

				if (children[i] == child)
					return true;

			}
			
			return false;

		}

		bool is_parent(node_id parent) {

			if (parent < 0)
				return false;

			for(size_t i = 0; i < parents.size(); i++) {

				if (parents[i] == parent)
					return true;

			}
			
			return false;

		}

		bool delete_child(node_id child) {

			if(child < 0)
				return false;

			vector<node_id>::iterator it = children.begin();

			for( ; it != children.end(); it++) {

				if (*it == child)
					break;

			}

			if(it != children.end()) {


				children.erase(it);
				return true;

			} else {

				return false;

			}

		}

		bool delete_parent(node_id parent) {

			if(parent < 0)
				return false;

			vector<node_id>::iterator it = parents.begin();

			for ( ; it != parents.end(); it++) {

				if (*it == parent)
					break;

			}

			if (it != parents.end()) {
			
				parents.erase(it);
				return true;

			} else {

				return false;

			}

		}

		void clear_children() {

			children.clear();

		}

		void clear_parents() {

			parents.clear();

		}

		void clear() {

			sequence.clear();
			children.clear();
			parents.clear();

		}
		~Node() { }

	};


//template <typename T>
class SplicingGraph {

private:

	class SortedNode {
		
	public:

		SortedNode(vector<Node>& sn) : sorted_node(sn) {};

		bool operator() (const node_id i, const node_id j){

			return (sorted_node[i].sequence.length() < sorted_node[j].sequence.length() ||
				(sorted_node[i].sequence.length() == sorted_node[j].sequence.length() && i < j));

		}

	private:

		vector<Node>& sorted_node;

	};

public:

	vector<Node> node_set;
	size_t node_sum;
	vector<int> node_order;
        vector<node_id> forward_branches;
	//set<pair_t> Inhibit_edges;
	//vector<vector<int>> Pair_edges;

	SplicingGraph();
	//~SplicingGraph();
	size_t get_node_sum();
	void set_parents();
	int add_node(Node& node);
	bool build(ReadHash& read_hash, KmerHash& kmer_hash, int seed, vector<pair<string,float> >& transcripts);
	void refine_trunk(ReadHash& read_hash, KmerHash& kmer_hash);
	bool is_trunk(ReadHash& read_hash, const string& trunk);
	//bool dele_node(const node_id p);
	void refine_forward(ReadHash& read_hash, KmerHash& kmer_hash, node_id p, int& num);
	void refine_reverse(ReadHash& read_hash, KmerHash& kmer_hash, node_id p, int& num);
	void refine_graph(ReadHash& read_hash, KmerHash& kmer_hash);
	void map_reads_to_sequence(ReadHash& kmer_hash, const string& sequence, vector<pair<int, int> >& reads_pos);
	string forward_extend(ReadHash& read_hash, KmerHash& kmer_hash, const string& seed_contig, int color_tag = -1);
	string reverse_extend(ReadHash& read_hash, KmerHash& kmer_hash, const string& seed_contig, int color_tag = -1);	
	bool forward_extend(ReadHash& read_hash, KmerHash& kmer_hash, string& contig, const string& stop_seq, int color_tag = -1);
	bool reverse_extend(ReadHash& read_hash, KmerHash& kmer_hash, string& contig, const string& stop_seq, int color_tag = -1);	
	bool pair_support(ReadHash& read_hash, const string& sequence, node_id p);
	bool is_circle(node_id p, node_id q, set<node_id>& checked);
	void trim_graph(ReadHash& read_hash); 
	void set_coverage_of_nodes(ReadHash& read_hash);
	float compute_coverage(ReadHash& read_hash, const string& sequence);
	string get_edge_sequence(node_id s, node_id t);
	void get_coverage_of_edges(ReadHash& read_hash, map<pair_t, double>& edge_coverages, vector<pair_t>& edges);
	void get_coverage_of_edges(ReadHash& read_hash, map<pair_t, double>& edge_coverages);
	void topological_sort();
	void dfs_visit(node_id i, vector<int>& node_color);  	
	int get_total_amount_of_kmers();
	void recover_reads();
	void set_reads_tag(ReadHash& read_hash, const string& sequence, int tag);
	void set_reads_tag(ReadHash& read_hash, const string& sequence, vector<int>& map_reads, int tag);
	void get_transcripts(ReadHash& read_hash, vector<pair<string,float> >& transcripts);
	void get_transcripts1(ReadHash& read_hash, vector<pair<string,float> >& transcripts);
	void check_transcripts(ReadHash& read_hash, vector<pair<string,float> >& transcripts);
	//void get_transcripts_by_pair(ReadHash& read_hash, vector<pair<string,float> >& transcripts); 
	int get_read_id(node_id p, int pos);
	void get_edge_reads(ReadHash & read_hash, node_id source, node_id target, vector<int>& edge_reads);
	void solve_milp(vector<pair<node_id, float> >& out_edge, vector<pair<node_id, float> >& in_edge, vector<float>& adj, vector<vector<int> >& in_candi, vector<vector <float> >& in_priority);
	void set_reads_pos_in_node(ReadHash& read_hash, node_id p);
	void branch_extend_by_coverage(ReadHash& read_hash, KmerHash& kmer_hash);
	void forward_extend_by_coverage(ReadHash& read_hash, KmerHash& kmer_hash, node_id p);
	void reverse_extend_by_coverage(ReadHash& read_hash, KmerHash& kmer_hash, node_id p);
	void set_forward_check_pos(ReadHash& read_hash, node_id p);
	void set_reverse_check_pos(ReadHash& read_hash, node_id p);
	void refine_forward_trunk(ReadHash& read_hash, KmerHash& kmer_hash);
	void refine_reverse_trunk(ReadHash& read_hash, KmerHash& kmer_hash);
	void set_trunk_check_pos(ReadHash& read_hash,node_id p);
	bool refine_forward_by_pair(ReadHash& read_hash, KmerHash& kmer_hash, node_id p);
	bool refine_reverse_by_pair(ReadHash& read_hash, KmerHash& kmer_hash, node_id p);
	~SplicingGraph() { }
};



#endif
