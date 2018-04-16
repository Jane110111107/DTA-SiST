
#ifndef SPLICINGGRAPH_H
#define SPLICINGGRAPH_H

//extern "C" {
//#include "LingoProgram.h"
//}

#include"KmerUtility.h"
#include "ReadHash.h"
#include"SuffixTree.h"
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



int get_compatible_len(const string& str1, const string& str2);
bool is_seq_similar(const string& str1, const string& str2, char mode = 'F', double sim_error = 0.35);
bool is_aligned(const string& str1, const string& str2, char mode = 'F', int tag = 0);
bool compatible(const string& str1, const string& str2);




class Node {

	public:

		string sequence;
		vector<int> parents;
		vector<int> children;
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

		bool add_child(int child) {

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

		bool add_parent(int parent) {

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

		bool is_child(int child) {

			if (child < 0)
				return false;

			for(int i = 0; i < children.size(); i++) {

				if (children[i] == child)
					return true;

			}
			
			return false;

		}

		bool is_parent(int parent) {

			if (parent < 0)
				return false;

			for(size_t i = 0; i < parents.size(); i++) {

				if (parents[i] == parent)
					return true;

			}
			
			return false;

		}

		bool delete_child(int child) {

			if(child < 0)
				return false;

			vector<int>::iterator it = children.begin();

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

		bool delete_parent(int parent) {

			if(parent < 0)
				return false;

			vector<int>::iterator it = parents.begin();

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


class SplicingGraph {

public:

	vector<Node> node_set;
	size_t node_sum;
	vector<int> node_order;
        vector<int> forward_branches;

	SplicingGraph();

	size_t get_node_sum();
	void set_parents();
	int add_node(Node& node);
	bool build(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree, int seed, vector<pair<string,float> >& transcripts);
	//string forward_extend(SuffixTree& right_tree, const string& seed_contig);
	//string reverse_extend(SuffixTree& left_tree, const string& seed_contig);
	string forward_extend(SuffixTree& right_tree, int seed);
	string reverse_extend(SuffixTree& left_tree, int seed);
	void set_reads_tag(ReadHash& read_hash, const string& sequence, int tag);
	void set_reads_tag(ReadHash& read_hash, const string& sequence, vector<int>& map_reads, int tag);
	bool is_trunk(ReadHash& read_hash, const string& trunk);
	void refine_forward_trunk(ReadHash& read_hash, SuffixTree& right_tree);
	void refine_reverse_trunk(ReadHash& read_hash, SuffixTree& left_tree);
	void set_reads_pos_in_node(ReadHash& read_hash, int p);
	bool refine_reverse_by_pair(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree, int p);
	bool refine_forward_by_pair(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree, int p);
	//bool forward_extend(SuffixTree& right_tree, string& contig, const string& stop_seq);
	//bool reverse_extend(SuffixTree& left_tree, string& contig, const string& stop_seq);
	bool forward_extend(SuffixTree& right_tree, string& contig, const string& stop_seq,int seed);
	bool reverse_extend(SuffixTree& left_tree, string& contig, const string& stop_seq, int seed);
	void branch_extend_by_coverage(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree);
	void forward_extend_by_coverage(ReadHash& read_hash, SuffixTree& right_tree, int p);
	void reverse_extend_by_coverage(ReadHash& read_hash, SuffixTree& left_tree, int p);
	void set_forward_check_pos(ReadHash& read_hash, int p);
	void set_reverse_check_pos(ReadHash& read_hash, int p);	
	bool pair_support(ReadHash& read_hash, const string& sequence, int p);
	bool is_circle(int p, int q, set<int>& checked);
	void refine_graph(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree);
	void topological_sort();
	void dfs_visit(int i, vector<int>& node_color); 
	void set_coverage_of_nodes(ReadHash& read_hash);
	float compute_coverage(ReadHash& read_hash, const string& sequence);
	void trim_graph(ReadHash& read_hash); 
	string get_edge_sequence(int s, int t);
	void get_coverage_of_edges(ReadHash& read_hash, map<pair_t, double>& edge_coverages, vector<pair_t>& edges);
	void get_coverage_of_edges(ReadHash& read_hash, map<pair_t, double>& edge_coverages);
	void get_transcripts_by_all(ReadHash& read_hash, vector<pair<string,float> >& transcripts);
	void get_all_path(int node, vector<int>& path_node, vector<pair<string,float> >& transcripts);
	void get_transcripts_by_cov_and_len(ReadHash& read_hash, vector<pair<string,float> >& transcripts);
	void recover_reads();
	void check_transcripts(ReadHash& read_hash, vector<pair<string,float> >& transcripts);
	void map_reads_to_sequence(ReadHash& kmer_hash, const string& sequence, vector<pair<int, int> >& reads_pos);

	~SplicingGraph() { }
};



#endif
