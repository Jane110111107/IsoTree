// This file was modified from a file named common.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpack Software LICENSE.



#include "GeneralSet.h"
#include <cstring>
#include <cstdlib>
#include <errno.h>

using namespace std;


vector<int> data_tag(1000000);
vector<string> data(1000000);
int max_read_id = 0;
double g_max_mapped_error = 0.02;
int g_mode = 1;
int g_kmer_length = 25;
bool g_help = false;
bool g_debug = false;
bool g_is_paired_end = true; //false;
int g_fr_strand = 1;
bool g_double_stranded_mode = false; //false;
int g_min_same_len = 25;
int g_max_same_len = 0;
int g_max_check_len = 200;
int g_min_read_sum = 1;
size_t g_read_length = 75;
int g_tag_length = 5;
int g_same_length = 10;
double g_tolerance_value = 0.4;
float g_min_pair_ratio = 0.5;

// assemble
int g_refine_check_len = 30;
int g_min_kmer_coverage = 1;//1;
float g_min_kmer_entropy = 0.0f;
float g_min_trans_cov = 0.001;
int g_min_seed_coverage = 2;
float g_min_seed_entropy = 1.5f;
int g_min_junction_coverage = 2;
int g_min_average_coverage = 2;
int g_min_reads_span_junction = 2;
int g_min_reads_support_branch = 2;
int g_min_anchor_length = 21;
float g_min_ratio_non_error = 0.04f;
float g_min_ratio_branch = 0.1f;
float g_min_ratio_welds = 0.04f;
float g_min_ratio_in_out = 0.1f;
int g_min_exon_length = 80;
int g_min_trunk_length = 200;
int g_min_kmers_per_graph = 276;
int g_pair_gap_length = 200;
int g_max_pair_gap_length = 500;
string g_reads_file = "";
string g_left_file = "";
string g_right_file = "";
string out_dir = "./SplicingGraphs/"; //change by Jane
string sg_list = "splicing_graph_name.list"; //change by Jane
int g_interval = 20;

// path search
int CPU = 6;
string rg_file = "";
int g_min_transcript_length = 200;
string output_filename = "transcripts.fasta";

// path search
