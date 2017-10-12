# IsoTree
IsoTree: A New Framework for De Novo Transcriptome Assembly from RNA-seq Reads

********************Description************************

IsoTree is an efficient de novo trascriptome assembler for RNA-Seq data. It can assemble transcripts from RNA-Seq reads (in fasta format). Unlike most of de novo assembly methods that build de Bruijn graph or splicing graph by connecting k-mers which are sets of overlapping substrings generated from reads, IsoTree constructs splicing graph by connecting reads directly. For each splicing graph, IsoTree applies an iterative scheme of mixed integer linear program to build a prefix tree, called isoform tree. Each path from the root node of the isoform tree to a leaf node represents a plausible transcript candidate which will be pruned based on the information of pair-end reads. 

**********************Install**************************

$gcc -I /path/to/boost/include -o isotree GeneralSet.cpp KmerUtility.cpp ReadUtility.cpp ReadHash.cpp TreeStruct.cpp SplicingGraph.cpp IsoTree.cpp -pthread -lirc -lstdc++ -limf
$export PATH=/path/to/isotree:$PATH

***********************Uasge***************************

quik start:

$isotree -t 1 -2 100 -k 20 --left sim100bp_1.fa --right sim100bp_2.fa --max_same_len 99 --min_same_len 20 --fr_strand 1
(where sim100bp_1.fa and sim100bp_2.fa are simulated paired-end reads ())

detail:

~ Required: ~ 
	 -t <int>: type of reads: 1: single-end reads,  2: paired-end reads.
	 -l <int>: read length. 
If pair end reads: 
	 --left <string>: left reads file name (.fasta). 
	 --right <string>: right reads file name (.fasta). 
If single end reads: 
	 --singlefile <string>: reads file name (.fasta).
 ~ Options: ~
	-k <int>: length of kmer, default 25. 
	-o <string>: name of drectory for output
	-h : help information. 
	--min_same_len <int>: the minimum overlap length, default k. 
	--max_same_len <int>: the maximum overlap length, default read_length -1.
	--tolerance_value <float>: the value of epsilon, default 0.35. 
	--min_trans_len <int>: the minimum length of transcript, default 200.
	--min_exon_len <int>: the minimum length of exon, default 80.
  --mode <int>: method of constructing splicing graph. 1: both trunk extend and branch extend, 2: only trunk extend, default 1. 
  --double_stranded_mode: indicate the pair-end read is double stranded mode" << endl
	--fr_strand <int>: only used for pair-end reads. 1: --1--> <--2--  2: <--1-- --2-->  3: --1--> --2--> or <--1-- <--2--, default 1.
		
