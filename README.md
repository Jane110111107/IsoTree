# IsoTree
IsoTree: A New Framework for De Novo Transcriptome Assembly from RNA-seq Reads

version 1.1

** Description **

IsoTree is an efficient de novo trascriptome assembler for RNA-Seq data. It can assemble transcripts from RNA-Seq reads (in fasta format). Unlike most of de novo assembly methods that build de Bruijn graph or splicing graph by connecting k-mers which are sets of overlapping substrings generated from reads, IsoTree constructs splicing graph by connecting reads directly. For each splicing graph, IsoTree applies an iterative scheme of mixed integer linear program to build a prefix tree, called isoform tree. Each path from the root node of the isoform tree to a leaf node represents a plausible transcript candidate which will be pruned based on the information of pair-end reads. 

** Install ** 
$export PATH=/path/tp/boost/include:$PATH
$g++ -o isotree GeneralSet.cpp KmerUtility.cpp KmerHash.cpp ReadUtility.cpp ReadHash.cpp TreeStruct.cpp SplicingGraph.cpp IsoformTree.cpp
$export PATH=/path/to/isotree:$PATH

** Uasge **

quik start:

$isotree -t 2 -l 100 -k 20 --left reads_1.fa --right reads_2.fa --max_same_len 99 --min_same_len 20 --fr_strand 1
(where -l 100 represents the length of reads and -k 20 corresponds to the length of k-mer)

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
	-h : help information. 
	--min_same_len <int>: the minimum overlap length, default k. 
	--max_same_len <int>: the maximum overlap length, default read_length -1.
	--tolerance_value <float>: the value of epsilon, default 0.35. 
	--min_trans_len <int>: the minimum length of transcript, default 200.
	--min_exon_len <int>: the minimum length of exon, default 80.
  --mode <int>: method of constructing splicing graph. 2: both trunk extend and branch extend, 1: only trunk extend, default 2. 
  --double_stranded_mode: indicate the pair-end read is double stranded mode" << endl
	--fr_strand <int>: only used for pair-end reads. 1: --1--> <--2--  2: <--1-- --2-->  3: --1--> --2--> or <--1-- <--2--, default 1.


** Output **

The IsoTree will output a fasta file named transcriptome.fa that contains all the possible transcripts assembled by IsoTree.


** Simulated data **

We generated total of 21 datasets of simulated pair-end reads from 100 isoform transcripts（stored in ref_trans.fa and ref_trans.gtf）originated from 41 different genes in chromosome 1 (CRCh38.83, NCBI) using FluxSimulator. Each simulated dataset contains 0.1 million paired-end reads. For example, sim100bp_1.fa and sim100bp_2.fa contain 0.1 million paired-end reads with length 100bp. 
Download link http://pan.baidu.com/s/1misseHa.




** Contact information **

If you have any questions or concerns, please send email to zhaojin_cc@163.com.
	

