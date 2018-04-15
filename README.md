# DTAST
DTAST: a novel radical framework for de novo transcriptome assembly based on suffix trees

** Description **
The RNA-seq has revolutionized our ability to transcripts
reconstruction. A growing number of strategies have been developed to
solve the transcriptome assembly problems, but their performances still
need to be improved. In this article, we develop a novel radical framework
for de novo transcriptome assembly based on suffix trees, called
DTAST. DTAST first extends contigs by reads that keep the longest
overlaps with the contigs’ terminuses. These reads can be found in the
linear time of the longest overlaps length through a well-designed suffix
tree structure. Then, DTAST constructs splicing graphs based on contigs
for each gene locus. Finally, DTAST proposes two strategies to extract
transcript-representing paths: depth-first enumeration strategy and
the hybrid strategy based on length and coverage. We implemented the
above two strategies and compared them with the state-of-the-art de novo
assemblers on both simulated and real datasets. Experimental results
showed that DTAST performs more competitive than the other compared
de novo assemblers especially when the read sequence is long.

** Uasge **
quik start:
$suffixtree -t 2 -l 100 -k 20 --left reads_1.fa --right reads_2.fa --max_same_len 99 --min_same_len 20 --fr_strand 1 (where -l 100 represents the length of reads and -k 20 corresponds to the length of k-mer)
detail:
~ Required: ~ -t : type of reads: 1: single-end reads, 2: paired-end reads. -l : read length. If pair end reads: --left : left reads file name (.fasta). --right : right reads file name (.fasta). If single end reads: --singlefile : reads file name (.fasta). ~ Options: ~ -k : length of kmer, default 25. -h : help information. --min_same_len : the minimum overlap length, default k. --max_same_len : the maximum overlap length, default read_length -1. --tolerance_value : the value of epsilon, default 0.35. --min_trans_len : the minimum length of transcript, default 200. --min_exon_len : the minimum length of exon, default 80. --mode : strategy to extract transcript-representing path. 1: depth-first enumeration strategy 2: the hybrid strategy based on length and coverag, default 1. --double_stranded_mode: indicate the pair-end read is double stranded mode" << endl --fr_strand : only used for pair-end reads. 1: --1--> <--2-- 2: <--1-- --2--> 3: --1--> --2--> or <--1-- <--2--, default 1.

** Output **
DTAST will output a fasta file named transcriptome.fa that contains all the possible transcripts assembled by DTAST.

** Contact information **
If you have any questions or concerns, please send email to zhaojin_cc@163.com.
