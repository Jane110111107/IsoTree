// This file was modified from a file named common.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpack Software LICENSE.



#ifndef GENERALSET_H
#define GENERALSET_H



#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>


using namespace std;

extern vector<int> data_tag;
extern vector<string> data;
extern int max_read_id;
extern double g_max_mapped_error;
extern int g_mode;
extern int g_kmer_length;
extern bool g_help;
extern bool g_debug;
extern bool g_double_stranded_mode;
extern bool g_is_paired_end;
extern string out_dir;
extern int g_fr_strand;
extern int g_min_same_len;
extern int g_max_same_len;
extern int g_max_check_len;

extern int g_min_read_sum;
extern int g_refine_check_len;
extern size_t g_read_length;
extern int g_tag_length;
extern int g_same_length;
extern double g_tolerance_value;
extern float g_min_pair_ratio;
// assemble options
extern string g_reads_file; //
extern string g_left_file; //
extern string g_right_file; //
extern int g_min_kmer_coverage;
extern float g_min_kmer_entropy;
extern float g_min_trans_cov;
extern int g_min_seed_coverage;
extern float g_min_seed_entropy;
extern int g_min_junction_coverage;
extern int g_min_average_coverage;
extern int g_min_anchor_length;
extern int g_min_reads_span_junction;
extern int g_min_reads_support_branch;
extern float g_min_ratio_non_error;
extern float g_min_ratio_branch;
extern float g_min_ratio_welds;
extern float g_min_ratio_in_out;
extern int g_min_exon_length;
extern int g_min_trunk_length;
extern int g_pair_gap_length;
extern int g_max_pair_gap_length;
extern int g_min_kmers_per_graph;
extern bool g_double_stranded_mode;
extern string sg_list;
extern int g_interval;///////delete

// path search options
extern int CPU;
extern string rg_file;
extern int g_min_transcript_length;
extern string output_filename;



#define OPT_MIN_SAME_LEN		301
#define OPT_MAX_SAME_LEN		302
#define OPT_DOUBLE_STRANDED_MODE		303
#define OPT_FR_STRAND		304
#define	OPT_TOLERANCE_VALUE	305
#define OPT_MIN_TRANS_LEN		306
#define OPT_MIN_EXON_LEN		307
#define OPT_LEFT		308
#define OPT_RIGHT	309
#define OPT_SINGLEFILE			311
#define OPT_MODE	312



#endif
