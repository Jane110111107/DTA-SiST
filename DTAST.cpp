

#include "GeneralSet.h"
#include "ReadUtility.h"
#include "ReadHash.h"
#include "SuffixTree.h"
#include "SplicingGraph.h"
#include <iostream>
#include <ctime>
#include <errno.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <getopt.h>
#include <stdio.h>

using namespace std;


struct option opts[] = {
	{"read_length",   required_argument,   0,   'l'},
	{"out_dir",        required_argument,   0,   'o'},
	{"is_paired",     required_argument,         0,   't'},
	{"help",          no_argument,         0,   'h'},
	{"min_same_len",  required_argument,   0,   'k'},
	{"model",          required_argument,   0,   OPT_MODEL},
	{"double_stranded_mode", no_argument,  0,   OPT_DOUBLE_STRANDED_MODE},
	{"fr_strand",     required_argument,   0,   OPT_FR_STRAND},
	{"tolerance_value",required_argument,  0,   OPT_TOLERANCE_VALUE},
	{"min_trans_len", required_argument,   0,   OPT_MIN_TRANS_LEN},
	{"min_exon_len",  required_argument,   0,   OPT_MIN_EXON_LEN},
	{"left",          required_argument,   0,   OPT_LEFT},
	{"right",         required_argument,   0,   OPT_RIGHT},
	{"singlefile",    required_argument,   0,   OPT_SINGLEFILE},
	{0,0,0,0}

};


string usage() {

	stringstream usage_info;
	usage_info
		<< endl
		<< "===============================================================================" << endl
		<< " IsoTree Usage " << endl
		<< "===============================================================================" << endl
		<< " ** Required: **" << endl
		<< "  -t <int>: type of reads: 1: single-end reads,  2: paired-end reads." << endl
		<< " If pair end reads: " << endl
		<< "  --left <string>: left reads file name (.fasta). " << endl
		<< "  --right <string>: right reads file name (.fasta). " << endl
		<< " If single end reads: " << endl
		<< "  -singlefile <string>: reads file name (.fasta). " << endl
		<< " ** Options: **" <<endl
		<< "  -l <int>: read length" << endl
		<< "  -k <int>: min overlap length of two reads, default 25. " << endl
		<< "  -o <string>: name of drectory for output, default: ./Isotree_Out/ " << endl
		<< "  -h : help information. " << endl
		<< "  --tolerance_value <float>: the value of epsilon, default 0.35" << endl
		<< "  --min_trans_len <int>: the minimum length of transcript, default 200." << endl
		<< "  --min_exon_len <int>: the minimum length of exon, default 80. " << endl
        << "  --mode <int>: method of constructing splicing graph. 1: both trunk extend and branch extend, 2: only trunk extend, default 1. " << endl
		<< "  --double_stranded_mode: indicate the pair-end read is double stranded mode" << endl
		<< "  --fr_strand <int>: only used for pair-end reads. 1: --1--> <--2--  2: <--1-- --2-->  3: --1--> --2--> or <--1-- <--2--, default 1. " << endl
		<< "===============================================================================" << endl
		<< endl;

	return usage_info.str();

}


int parse_options(int argc, char* argv[]) {

	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, "k:l:o:t:h", opts, &option_index);
		switch (next_option) {
		case -1:
			break;
		case 'k':
			g_min_same_len = atoi(optarg);
			break;
		case 'l':
			g_read_length = atoi(optarg);
			break;
		case 'o':
			out_dir = optarg;
			break;
		case 't':
			if (atoi(optarg) == 1)			
				g_is_paired_end = false;
			break;
		case 'h':
			g_help = true;
			break;
		case OPT_MODEL:
			g_mode = atoi(optarg);
			break;
		case OPT_DOUBLE_STRANDED_MODE:
			g_double_stranded_mode = true;
			break;
		case OPT_FR_STRAND:
			g_fr_strand = atoi(optarg);
			break;
		case OPT_TOLERANCE_VALUE:
			g_tolerance_value = atof(optarg);
			break;
		case OPT_MIN_TRANS_LEN:
			g_min_transcript_length = atoi(optarg);
			break;
		case OPT_MIN_EXON_LEN:
			g_min_exon_length = atoi(optarg);
			break;
		case OPT_LEFT:
			g_left_file = optarg;
			break;
		case OPT_RIGHT:
			g_right_file = optarg;
			break;
		case OPT_SINGLEFILE:
			g_reads_file = optarg;
                        break;
		default:
			cout << usage();
			exit(1);
		}

	} while (next_option != -1);

	if (g_help) {
		cout << usage();
		exit (1);
	}

	if (g_is_paired_end && (g_right_file == "" || g_left_file == "")) {
		cout << "Error! Please check --left and --right options." << endl;
		exit(1);
	}

	if (g_read_length == 0) {
		cout << "Error! Please add the read length by -l option." << endl;
		exit(1);
	}


	if (!g_is_paired_end && g_reads_file == "") {
		cout << "Error! Please check --singlefile option." << endl;
		exit(1);
	}

	if (g_min_same_len > 32) {
		cout << "Error: the minimum overlap length should shorter than 32." << endl;
		exit(1);
	}


	if (g_fr_strand != 1 && g_fr_strand != 2 && g_fr_strand != 3) {
		cout << "Error: --fr_strand can only be 1, 2, or 3" << endl;
		exit(1);
	}

	return 0;

}


void assemble(){
	
	ReadHash read_hash;
	read_hash.get_read_hash();
	SuffixTree left_tree;
	left_tree.get_left_tree();
	SuffixTree right_tree;
	//right_tree.get_right_tree();
	left_tree.get_right_tree(right_tree);

	const string& transcriptome_name = "transcriptome.fa";
	fstream transcriptome_file;
	transcriptome_file.open(transcriptome_name.c_str(), fstream::out);
	int gene_id = 0;
	if (!transcriptome_file.is_open()) {
		cout <<"File " << transcriptome_name.c_str() <<" can't be opened. " << endl;
		exit(1);
	}

	time_t s_time = time(NULL);
	const int data_size = data.size();
	size_t splicing_graph_id = 0;
	cout << "Assembling transcripts..." << endl;
	cout <<"Done 0 : " << data_size << endl;

	for (unsigned int i = 0; i < data_size; ++i) {

		if (data_tag[i] != -5) 
			continue;

		SplicingGraph splicing_graph;			
		//cout << "Building splicing graph" << splicing_graph_id << " ..." << endl;
		vector<pair<string,float> > transcripts;		
		if (splicing_graph.build(read_hash, left_tree, right_tree, i, transcripts)) {
			int sg_num = 0;
			if (transcripts.size() == 0)
				continue;
			for (size_t k =0; k < transcripts.size(); k++){				
				sg_num++;
				gene_id++;
				transcriptome_file<<">trans" << gene_id << "_sg" << splicing_graph_id << "_" << sg_num << "  len = " << transcripts[k].first.length() << "  cov = " << transcripts[k].second << "  sequence =" << endl;
				transcriptome_file << transcripts[k].first <<endl;
			}				
			if (sg_num > 0) {
				//cout << "Building splicing graph " << splicing_graph_id << " succeed!" << " total node sum: " << splicing_graph.get_node_sum() << endl;
				splicing_graph_id++;
			}	

		} 

		printf("\033[1ADone \%d : %d \n", i, data_size);
//cout << "Done " << i << ":" << data_size << endl;
		
	}
	printf("\033[1ADone \%d : %d \n", data_size, data_size);
	transcriptome_file.close();
	time_t e_time = time(NULL);		
	cout << "Total assembled " << gene_id << " transcripts (elapsed time: "  << (e_time - s_time) << " s)" << endl;


}


int main(int argc, char* argv[]){

	vector<string> input_data;
	int parse_ret = parse_options(argc,argv);
	if (parse_ret)
		return parse_ret;

	time_t s_time = time(NULL);
	if (g_is_paired_end) {
		//const int file_size = count_reads_file(g_left_file);
		//data.reserve(2*file_size+1);
//cout << data_tag.capacity() << endl;
		//data_tag.reserve(2*file_size+1);
//cout << data_tag.capacity() << endl;
		if (g_double_stranded_mode) {								
			load_reads(g_left_file,input_data, false);								
			load_reads(g_right_file,input_data,false);               
		} else {              
			if (g_fr_strand == 1) {//--1-->  <--2--				
				load_reads(g_left_file, input_data, false);				
				load_reads(g_right_file, input_data, true);
            }                	
			if (g_fr_strand == 2) {//<--1-- --2-->				
				load_reads(g_left_file, input_data, true);				
				load_reads(g_right_file,input_data, false);               	 
			}
			if (g_fr_strand == 3) {//--1--> --2-->  or --2--> --1-->				
				load_reads(g_left_file, input_data, false);				
				load_reads(g_right_file, input_data, false);              	
			}               
		}

	} else {
		const int file_size = count_reads_file(g_reads_file);
		//data.reserve(file_size+1);
		//data_tag.reserve(file_size+1);
		load_reads(g_reads_file,input_data, false);

	}
	time_t e_time = time(NULL);
	cout <<"Total load " << input_data.size() << " reads (elapsed time: " << (e_time - s_time) << " s)" << endl;

	g_read_length = input_data[0].length();
	if (g_max_same_len == 0) {
		g_max_same_len = g_read_length - 1;
	}

	if (g_min_same_len > g_read_length -1 ) {
		g_min_same_len = g_read_length - 1;
	}

	if (g_max_same_len > g_read_length -1) {
		g_max_same_len = g_read_length - 1;
	}

	delete_error_reads(input_data);
	assemble();
	cout << "Success! " << endl;
	return 1;

}

