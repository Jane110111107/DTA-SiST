#include "SuffixTree.h"
#include<iostream>

using namespace std;


trie_node* create_trie_node() {

	trie_node* pNode = new trie_node();
	for (int i = 0; i < 4; i++) {

		pNode->children[i] = NULL;
	}

	return pNode;

}


bool delete_trie_node(trie_node* p) {

	bool is_delete = true;
	if ((p->read_set).size() > 0) 
		is_delete = false;
	
	for (int i = 0; i < 4; ++i) {
		if (p->children[i] != NULL) {
			if (delete_trie_node(p->children[i])) {
				p->children[i] = NULL;
			}else {
				is_delete = false;
			}
		}
	}

	if (is_delete)
		delete p;

	return is_delete;	

}




bool SuffixTree::exists(kmer_int_type kmer) {

	suffix_tree_type_iterator it = find_kmer(kmer);
	return (it != suffix_tree.end());

}




void trie_insert(trie_node* root, const string& sequence, int r_d){
	//tag == 0: forward tag == 1: reverse
	
	trie_node* index = root;
	if (sequence.length() == 0) {
		(index -> read_set).push_back(r_d);
		return;
	} 	
	for (size_t j = 0; j < sequence.length(); ++j) {
		int pos = base_to_int(sequence[j]);
		index = index->children[pos];
		if (index == NULL)
			return;
	}

	(index->read_set).push_back(r_d);

}



void trie_insert_reverse(trie_node* root, const string& sequence, int r_d){
	//tag == 0: forward tag == 1: reverse
	
	trie_node* index = root;
	if (sequence.length() == 0) {
		(index -> read_set).push_back(r_d);
		return;
	} 	
	for (int j = sequence.length()-1; j >= 0; --j) {
		int pos = base_to_int(sequence[j]);
		index = index->children[pos];
		if (index == NULL)
			return;
	}

	(index->read_set).push_back(r_d);

}




trie_node* & SuffixTree::operator[](kmer_int_type kmer) {

	return suffix_tree[kmer];

}




void SuffixTree::get_right_tree() {

	const int data_size = data.size();

	cout << "Constructing right extension suffix tree..." << endl;

	time_t start_time = time(NULL);
	
	if (data.empty())
		return;
	if (g_read_length < g_min_same_len+1) {
		cout << "Error: the min overlap length is longer than the read length." << endl;
		return;
	}

	kmer_int_type kmer;
	trie_node* index;
	int pos = 0;

	for (size_t i = 0; i < data_size; ++i) {
		//if (data[i].length() < g_min_same_len + 1)
			//continue;
		const string& same_str = data[i].substr(g_read_length-g_min_same_len);
		kmer = kmer_to_int(same_str);
		if (!exists(kmer)) {
			index = create_trie_node();
			suffix_tree[kmer] = index;
		} else {
			index = suffix_tree[kmer];
		}
		const string& insert_str = data[i].substr(1, g_read_length-g_min_same_len-1);
		for (int j = insert_str.length()-1; j >= 0; --j) {
			pos = base_to_int(insert_str[j]);
			if (index->children[pos] == NULL) {
				index->children[pos] = create_trie_node();
			}
			index = index->children[pos];
//cout << j << " "<< insert_str << endl;
		}
	}
	
cout << "dsdfd" << endl;
	for (size_t i = 0; i < data_size; ++i) {
		for (int j = g_max_same_len; j >= g_min_same_len; --j) {
			const string& sequence = data[i].substr(0,j);
			const string& kmer = sequence.substr(sequence.length()-g_min_same_len);
			kmer_int_type kmer_int = kmer_to_int(kmer);
			if (!exists(kmer_int))
				continue;
			index = suffix_tree[kmer_int];
			const string& insert_str = sequence.substr(0, sequence.length()-g_min_same_len);
			trie_insert_reverse(index, insert_str, i);
		}	
	}

cout << "eefgf" << endl;

	for (size_t i = 0; i < data_size; ++i) {

		const string& same_str = data[i].substr(g_read_length-g_min_same_len);
		kmer = kmer_to_int(same_str);
		if (!exists(kmer)) 
			continue;
		index = suffix_tree[kmer];
		for (int j = 0; j < 4; ++j) {
			if (index->children[j] != NULL) {
				if (delete_trie_node(index->children[j])) {
					index->children[j] = NULL;
				}
			}
		}

		bool is_delete = true;
		if ((index->read_set).size() > 0)
			is_delete = false;
		else {
			for (int j = 0; j < 4; ++j) {
				if (index->children[j] != NULL)
					is_delete = false;
			}	
		}
		if (is_delete) {
			delete index;
			suffix_tree_type_iterator it = find_kmer(kmer);
			if (it != suffix_tree.end()) {
				it -> second = NULL;
     				suffix_tree.erase(it);
			}
		}
		
	}

	time_t end_time = time(NULL);
	cout << "The suffix tree for right extension has been constructed, total (elapsed time: " << (end_time - start_time) << " s)" << endl;

}



void SuffixTree::get_left_tree() {

	const int data_size = data.size();

	cout << "Constructing left extension suffix tree..." << endl;

	time_t start_time = time(NULL);
	
	if (data.empty())
		return;
	if (g_read_length < g_min_same_len+1) {
		cout << "Error: the min overlap length is longer than the read length." << endl;
		return;
	}

	kmer_int_type kmer;
	trie_node* index;
	int pos = 0;

	for (size_t i = 0; i < data_size; ++i) {

		const string& same_str = data[i].substr(0, g_min_same_len);
		kmer = kmer_to_int(same_str);
		if (!exists(kmer)) {
			index = create_trie_node();
			suffix_tree[kmer] = index;
		} else {
			index = suffix_tree[kmer];
		}
		const string& insert_str = data[i].substr(g_min_same_len, g_read_length-g_min_same_len-1);
		for (size_t j = 0; j < insert_str.length(); ++j) {
			pos = base_to_int(insert_str[j]);
			if (index->children[pos] == NULL) {
				index->children[pos] = create_trie_node();
			}
			index = index->children[pos];
		}
	}
	
	for (size_t i = 0; i < data_size; ++i) {
		for (size_t j = 1; j <= g_read_length-g_min_same_len; ++j) {
			const string& sequence = data[i].substr(j);
			const string& kmer = sequence.substr(0, g_min_same_len);
			kmer_int_type kmer_int = kmer_to_int(kmer);
			if (!exists(kmer_int))
				continue;
			index = suffix_tree[kmer_int];
			const string& insert_str = sequence.substr(g_min_same_len);
			trie_insert(index, insert_str, i);
		}	
	}

	for (size_t i = 0; i < data_size; ++i) {

		const string& same_str = data[i].substr(0, g_min_same_len);
		kmer = kmer_to_int(same_str);
		if (!exists(kmer)) 
			continue;
		index = suffix_tree[kmer];
		for (int j = 0; j < 4; ++j) {
			if (index->children[j] != NULL) {
				if (delete_trie_node(index->children[j])) {
					index->children[j] = NULL;
				}
			}
		}
//*
		bool is_delete = true;
		if ((index->read_set).size() > 0)
			is_delete = false;
		else {
			for (int j = 0; j < 4; ++j) {
				if (index->children[j] != NULL)
					is_delete = false;
			}	
		}
		if (is_delete) {
			delete index;
			suffix_tree_type_iterator it = find_kmer(kmer);
			if (it != suffix_tree.end()) {
				it -> second = NULL;
     				suffix_tree.erase(it);
			}
		}
//*/		
	}

	time_t end_time = time(NULL);
	cout << "The suffix tree for left extension has been constructed, total (elapsed time: " << (end_time - start_time) << " s)" << endl;

}


//*
void SuffixTree::get_right_tree(SuffixTree& right_tree) {

	cout << "Constructing right extension suffix tree..." << endl;
	time_t start_time = time(NULL);
	
	kmer_int_type kmer_int;
	const int data_size = data.size();
	trie_node* index =NULL;
	int pos = 0;

	for (int i = 0; i < data_size; ++i){
		const string& read = data[i];
		const string& kmer = read.substr(0,g_min_same_len);
		kmer_int = kmer_to_int(kmer);
		if (!exists(kmer_int))
			continue;
		index = suffix_tree[kmer_int];
		if (index == NULL)
			continue;
		if ((index->read_set).size() > 0) {
			right_tree.add_node_to_right_tree(kmer, i);
		}
		for (int j = g_min_same_len; j < g_max_same_len; ++j) {
			pos = base_to_int(read[j]);
			index = index -> children[pos];
			if (index == NULL)
				break;

			if ((index->read_set).size() > 0){
				right_tree.add_node_to_right_tree(read.substr(0, j+1),i);
			}
		}
	}


	time_t end_time = time(NULL);
	cout << "The suffix tree for right extension has been constructed, total (elapsed time: " << (end_time - start_time) << " s)" << endl;

}

//*/

void SuffixTree::add_node_to_right_tree(const string& sequence, int id) {

	if (sequence.length() < g_min_same_len) {
		return;
	}

	const string& same_str = sequence.substr(sequence.length()-g_min_same_len);
	kmer_int_type kmer = kmer_to_int(same_str);
	trie_node* index;
	if (!exists(kmer)) {
		index = create_trie_node();
		suffix_tree[kmer] = index;
	} else {
		index = suffix_tree[kmer];
	}
	
	if (sequence.length() == g_min_same_len) {
		(index->read_set).push_back(id);
		return;
	}

	for (int i = sequence.length()-g_min_same_len-1; i >= 0; --i) {
		int pos = base_to_int(sequence[i]);
		if (index->children[pos] == NULL) {
			index->children[pos] = create_trie_node();
		}
		index = index->children[pos];
	}
	(index->read_set).push_back(id);
}



size_t SuffixTree::get_node_sum() {

	trie_node* index;
	suffix_tree_type_iterator it;
	size_t total = 0;
	for (it = suffix_tree.begin(); it != suffix_tree.end(); ++it) {
		index = it -> second;
		dfs_search(total, index);
	}
	return total;

}


void SuffixTree::dfs_search(size_t& total, trie_node* index) {

	if (index == NULL)
		return;
	total = total + 1;
	for (int i = 0; i < 4; ++i) {
		dfs_search(total, index->children[i]);
	}

}










