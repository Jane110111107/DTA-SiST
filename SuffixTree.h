// This file was modified from a file named splicing_graph.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpacker Software LICENSE.



#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

#include "KmerUtility.h"
#include "GeneralSet.h"
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

#include <boost/unordered_map.hpp>


using namespace std;

typedef pair<int,int> pair_t;

struct trie_node {

	vector<int> read_set;
	trie_node *children[4];

};


trie_node* create_trie_node();
bool delete_trie_node(trie_node* p);
void trie_insert(trie_node* root, const string& sequence, int r_d, int tag);



class SuffixTree {

private:

	typedef  boost::unordered_map<kmer_int_type, trie_node* > suffix_tree_type;
	typedef  boost::unordered_map<kmer_int_type, trie_node* >::iterator suffix_tree_type_iterator;
	typedef  boost::unordered_map<kmer_int_type, trie_node* >::const_iterator suffix_tree_type_const_iterator;
	suffix_tree_type suffix_tree;

public:

	SuffixTree () { }
	trie_node* & operator[](kmer_int_type kmer);

	suffix_tree_type_iterator find_kmer(kmer_int_type kmer){
		return suffix_tree.find(kmer);
	}
	bool exists(const kmer_int_type kmer); //weather find the kmer in hash table
	void get_left_tree();
	void get_right_tree();
	void get_right_tree(SuffixTree& right_tree);
	void add_node_to_right_tree(const string& sequence, int id);
	size_t get_node_sum();
	void dfs_search(size_t& total, trie_node* index);	
};









#endif
