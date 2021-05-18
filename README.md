
# DTA-SiST

DTA-SiST: a novel radical framework for de novo transcriptome assembly based on suffix trees

** Description **

The RNA-seq has revolutionized our ability to transcripts reconstruction. A growing number of strategies have been developed to
solve the transcriptome assembly problems, but their performances still need to be improved. In this article, we develop a novel radical framework for de novo transcriptome assembly based on suffix trees, called DTAST. DTAST first extends contigs by reads that keep the longest overlaps with the contigs’ terminuses. These reads can be found in the linear time of the longest overlaps length through a well-designed suffix tree structure. Then, DTAST constructs splicing graphs based on contigs for each gene locus. Finally, DTAST proposes two strategies to extract transcript-representing paths: depth-first enumeration strategy and the hybrid strategy based on length and coverage. We implemented the above two strategies and compared them with the state-of-the-art de novo assemblers on both simulated and real datasets. Experimental results showed that DTA-SiST performs more competitive than the other compared de novo assemblers especially when the read sequence is long.



** Install ** 

$export PATH=/path/to/boost/include:$PATH 

$g++ -o suffixtree GeneralSet.cpp ReadUtility.cpp KmerUtility.cpp SuffixTree.cpp ReadHash.cpp SplicingGraph.cpp DTAST.cpp

$export PATH=/path/to/suffixtree:$PATH



** Uasge **

quik start:

$suffixtree -t 2 -l 100 -k 20 --left reads_1.fa --right reads_2.fa --max_same_len 99 --min_same_len 20 --fr_strand 1 (where -l 100 represents the length of reads and -k 20 corresponds to the length of k-mer)

detail:

~ Required: ~ -t : type of reads: 1: single-end reads, 2: paired-end reads. -l : read length. If pair end reads: --left : left reads file name (.fasta). --right : right reads file name (.fasta). If single end reads: --singlefile : reads file name (.fasta). ~ Options: ~ -k : length of kmer, default 25. -h : help information. --min_same_len : the minimum overlap length, default k. --max_same_len : the maximum overlap length, default read_length -1. --tolerance_value : the value of epsilon, default 0.35. --min_trans_len : the minimum length of transcript, default 200. --min_exon_len : the minimum length of exon, default 80. --mode : strategy to extract transcript-representing path. 1: depth-first enumeration strategy 2: the hybrid strategy based on length and coverag, default 1. --double_stranded_mode: indicate the pair-end read is double stranded mode" << endl --fr_strand : only used for pair-end reads. 1: --1--> <--2-- 2: <--1-- --2--> 3: --1--> --2--> or <--1-- <--2--, default 1.

** Output **

DTAST will output a fasta file named transcriptome.fa that contains all the possible transcripts assembled by DTAST.



** Simulated data **

We used FluxSimulator to simulate five samples with read lengths of 50bp, 75bp, 100bp, 125bp, and 150bp, respectively. The only difference between these five samples is the length of reads. Each sample contains 0.1 million paired-end reads that are generated from 100 isoform transcripts originated from 41 different genes in chromosome 1 (CRCh38.83, NCBI).

sim50bp_1.fa and sim50bp_2.fa contain total of 0.1 million paired-end reads with length of 50bp. 

sim75bp_1.fa and sim75bp_2.fa contain total of 0.1 million paired-end reads with length of 75bp. 

sim100bp_1.fa and sim100bp_2.fa contain total of 0.1 million paired-end reads with length of 100bp. 

sim125bp_1.fa and sim125bp_2.fa contain total of 0.1 million paired-end reads with length of 125bp. 

sim150bp_1.fa and sim150bp_2.fa contain total of 0.1 million paired-end reads with length of 150bp. 

reference.fa contains total of 100 reference transcripts.



** Contact information **

If you have any questions or concerns, please send email to zhaojin@mail.sdu.edu.cn.
