
#include "SplicingGraph.h"
#include<iostream>

using namespace std;


int get_compatible_len(const string& str1, const string& str2) {

	int compatible_len = 0;
	for (int i = g_min_same_len-1; i > 5; --i) {
		const string& comp1 = str1.substr(0, i);
		const string& comp2 = str2.substr(str2.length()-i);
		if (comp1 == comp2) {
			compatible_len = i;
			break;
		}
	}
	
	return compatible_len;

}


bool is_seq_similar(const string& str1, const string& str2, char mode, double sim_error) {

	int mismatch = 0;
	int length = str1.length() > str2.length() ? str2.length(): str1.length();

	if (length == 0) {

		if (str1.empty() && str2.empty()) 
			return true;
		else 
			return false;
		
	}

	if (mode == 'F') {

		for (int i = 0; i < length; ++i) {

			if (str1[i] != str2[i])
				mismatch++;

		}

	} else {

		for (int i = 0; i < length; ++i) {

			if (str1[str1.length()-i-1] != str2[str2.length()-i-1]) 
				mismatch++;

		}

	}

	if ((float)mismatch/length < sim_error) 
		return true;
	else
		return false;

}


bool is_aligned(const string& str1, const string& str2, char mode, int tag) {

	int mismatch = 0;
	int length = str1.length() > str2.length()? str2.length(): str1.length();

	if (mode == 'F') {	
		for(int i = 0; i < length; ++i) {	
			if(str1[i] != str2[i])						
				mismatch++;
		}
	} else {
		for (int i = 0; i < length; ++i) {
			if (str1[str1.length()-i-1] != str2[str2.length()-i-1])
				mismatch++;
		}
	}
	return (mismatch <= tag);

}


bool compatible(const string& str1, const string& str2) {

	for(unsigned int i = 0; i <= str2.length()-g_min_same_len; ++i) {

		const string& kmer = str2.substr(i,g_min_same_len);
		string::size_type start = str1.find(kmer);

		if(start != string::npos) {

			if (start > i)
				return is_aligned(str1.substr(start-i),str2,'F',2);
			else
				return is_aligned(str1,str2.substr(i-start),'F',2);

		}

	}
	
	return false;

}


SplicingGraph::SplicingGraph() {

	node_sum = 0;

}


size_t SplicingGraph::get_node_sum() {

	return node_sum;

}


void  SplicingGraph::set_parents() {

	for (int i = 0; i < node_sum; ++i) {

		if (!node_set[i].parents.empty())
			node_set[i].parents.clear();

	}

	for (size_t i = 0; i < node_sum; ++i) {

		vector<int>::iterator it;

		for (it = node_set[i].children.begin(); it != node_set[i].children.end(); ++it)
			node_set[*it].add_parent(i);

	}

}



int SplicingGraph::add_node(Node& node) {

	node_set.push_back(node);
	return (node_sum++);

}



bool SplicingGraph::build(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree, int seed, vector<pair<string,float> >& transcripts) {

	data_tag[seed] = -1;
	//cout << " Get trunk..." << endl;
	const string& left = reverse_extend(left_tree, seed);
	const string& right = forward_extend(right_tree, seed);
	string trunk = left;
	if (right.length() > g_read_length)
		trunk = left + right.substr(g_read_length);

	if(static_cast<int>(trunk.length()) < g_min_trunk_length) {
		//cout << "The trunk is too short!" << endl;
		set_reads_tag(read_hash, trunk, -4);		
		return false;
	}
//*
	if (!is_trunk(read_hash, trunk)) {
		set_reads_tag(read_hash, trunk, -3);		
		return false;
	}
//*/


	Node node(trunk);
	int p = add_node(node);
	//refine_forward_trunk(read_hash, right_tree);
	//refine_reverse_trunk(read_hash, left_tree);

	if (g_is_paired_end) {
		//if (refine_reverse_by_pair(read_hash,left_tree, right_tree,p) || refine_forward_by_pair(read_hash,left_tree,right_tree,p))
			set_reads_tag(read_hash, node_set[p].sequence, p);
			
	} else {
		set_reads_tag(read_hash, node_set[p].sequence, p);
	}

	//cout << "branch check and extend..." << endl;
	branch_extend_by_coverage(read_hash, left_tree, right_tree);

	if (g_is_paired_end)
		refine_graph(read_hash, left_tree, right_tree);

	topological_sort();

	if (g_mode == 1) {
		get_transcripts_by_all(read_hash, transcripts);
	} else {
		get_transcripts_by_cov_and_len(read_hash, transcripts);
	}
		
	if (transcripts.size() == 0)
		return false;

	return true;

}



string SplicingGraph::forward_extend(SuffixTree& right_tree, int seed) {

//cout << "forward extend..." << endl;

	string contig = data[seed];
	kmer_int_type kmer_int;
	trie_node* index;
	vector<vector<int> > candi_set;
	vector<int> candi_len;
	int overlap_len, pos;
	vector<int> read_vec;
	int read_id, mate_id, candi_id;
	bool is_extend = true;
	int contig_cov = data_cov[seed];
	int candi_cov = contig_cov*10000;

	while (is_extend) {

		const string& kmer = contig.substr(contig.length()-g_min_same_len);
		kmer_int = kmer_to_int(kmer);
		index = right_tree[kmer_int];
		if (index == NULL)
			break;
		candi_set.clear();
		candi_len.clear();
		overlap_len = g_min_same_len;
		//while ((index != NULL)&&overlap_len < g_read_length) {
		while (index != NULL) {
			if ( (index->read_set).size() > 0){
				read_vec = index -> read_set;
				candi_set.push_back(read_vec);
				candi_len.push_back(overlap_len);
			}
			++overlap_len;
			pos = base_to_int(contig[contig.length()-overlap_len]);
			index = index -> children[pos];
		}		
	
		if (candi_len.size() == 0)
			break;

		is_extend = false;

		if (g_is_paired_end) {
			for (int i = candi_set.size()-1; i >= 0; --i) {
				read_vec = candi_set[i];
				for (int j = 0; j < read_vec.size(); ++j) {
					read_id = read_vec[j];
					if (data_tag[read_id] > -2)
						continue;
					for (int k = 0; k < data_pair[read_id].size(); ++k) {
						mate_id = data_pair[read_id][k];
						if (data_tag[mate_id] > -2) {
							contig = contig + data[read_id].substr(candi_len[i]);
							data_tag[read_id] = -1;
							contig_cov = data_cov[read_id];
							is_extend = true;
							break;
						}
					}
					if (is_extend)
						break;
				}
				if(is_extend)
					break;				
			}

		}//if (g_is_paired_end)
		
		if (!is_extend) {
			for (int i = candi_set.size()-1; i >= 0; --i) {
				read_vec = candi_set[i];
				candi_id = -1;
				for (int j = 0; j < read_vec.size(); ++j) {
					read_id = read_vec[j];
					if (data_tag[read_id] > -2)
						continue;
					if (j == 0 || candi_cov < abs(data_cov[read_id] - contig_cov) ) {
						candi_id = read_id;
						candi_cov = abs(data_cov[read_id] - contig_cov);
					}
				}
				if (candi_id >= 0){
					contig = contig + data[candi_id].substr(candi_len[i]);
					data_tag[candi_id] = -1;
					contig_cov = data_cov[candi_id];
					is_extend = true;
					break;
				}				
			}
			
		} //if (!is_extend)		
	}

	return contig;
}




string SplicingGraph::reverse_extend(SuffixTree& left_tree, int seed) {

//cout << "reverse extend..." << endl;
	
	string contig = data[seed];
	kmer_int_type kmer_int;
	vector<vector<int> > candi_set;
	vector<int> candi_len;
	int overlap_len, pos;
	vector<int> read_vec;
	int read_id, mate_id, candi_id;
	bool is_extend = true;
	trie_node* index = NULL;
	int contig_cov = data_cov[seed];
	int candi_cov = contig_cov*10000;

	while (is_extend) {

		const string& kmer = contig.substr(0, g_min_same_len);
		kmer_int = kmer_to_int(kmer);
		if (!left_tree.exists(kmer_int))
			break;
		index = left_tree[kmer_int];
		if (index == NULL)
			break;
		candi_set.clear();
		candi_len.clear();
		overlap_len = g_min_same_len;
		//while ((index != NULL)&&overlap_len<g_read_length) {
		while (index != NULL) {

			if ( (index->read_set).size() > 0){
				read_vec = index -> read_set;
				candi_set.push_back(read_vec);
				candi_len.push_back(overlap_len);
			}
			pos = base_to_int(contig[overlap_len]);
			index = index -> children[pos];
			++overlap_len;
		}	
	
		if (candi_len.size() == 0)
			break;

		is_extend = false;

		if (g_is_paired_end) {
			for (int i = candi_set.size()-1; i >= 0; --i) {
				read_vec = candi_set[i];
				for (int j = 0; j < read_vec.size(); ++j) {
					read_id = read_vec[j];
					if (data_tag[read_id] > -2)
						continue;
					for (int k = 0; k < data_pair[read_id].size(); ++k) {
						mate_id = data_pair[read_id][k];
						if (data_tag[mate_id] > -2) {
							contig = data[read_id].substr(0, g_read_length-candi_len[i]) + contig;
							data_tag[read_id] = -1;
							contig_cov = data_cov[read_id];
							is_extend = true;
							break;
						}
					}
					if (is_extend)
						break;
				}
				if(is_extend)
					break;				
			}

		}//if (g_is_paired_end)
		
		if (!is_extend) {
			for (int i = candi_set.size()-1; i >= 0; --i) {
				read_vec = candi_set[i];
				candi_id = -1;
				for (int j = 0; j < read_vec.size(); ++j) {
					read_id = read_vec[j];
					if (data_tag[read_id] > -2)
						continue;
					if (j == 0 || candi_cov < abs(data_cov[read_id] - contig_cov) ) {
						candi_id = read_id;
						candi_cov = abs(data_cov[read_id] - contig_cov);
					}
				}
				if (candi_id >= 0){
					contig = data[candi_id].substr(0, g_read_length-candi_len[i]) + contig;
					is_extend = true;
					data_tag[candi_id] = -1;
					contig_cov = data_cov[candi_id];
					break;
				}				
			}
		} //if (!is_extend)		
	}

	return contig;
}



void SplicingGraph::set_reads_tag(ReadHash& read_hash, const string& sequence, int tag) {

	if (sequence.length() < g_read_length)
		return;
	read_int_type read_int;
	int read_id;
	for (int i = 0; i <= sequence.length() - g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_id = read_hash.find_read(read_int);
		if (read_id < 0 || data_tag[read_id] == -2)
			continue;
		data_tag[read_id] = tag;
	}
}


void SplicingGraph::set_reads_tag(ReadHash& read_hash, const string& sequence, vector<int>& map_reads, int tag) {

	if (sequence.length() < g_read_length)
		return;

	read_int_type read_int;
	int read_id;
	for (int i = 0; i <= sequence.length() - g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_id = read_hash.find_read(read_int);
		if (read_id < 0 || data_tag[read_id] == -2)
			continue;
		data_tag[read_id] = tag;
		map_reads.push_back(read_id);
	}
}



bool SplicingGraph::is_trunk(ReadHash& read_hash, const string& trunk) {

	if (trunk.length() < g_read_length)
		return true;

	int read_id;
	int used_sum = 0;
	int total_sum = 0;
	for (int i = 0; i <= trunk.length()-g_read_length; ++i) {
		const string& read = trunk.substr(i, g_read_length);
		read_int_type read_int = get_read_int_type(read);
		read_id = read_hash.find_read(read_int);
		if (read_id < 0)
			continue;
		++total_sum;
		if (data_tag[read_id] == -3)
			++used_sum;		
	}
	if (used_sum > 0.8*total_sum && trunk.length() < 1000)
		return false;
	if (used_sum > 0.9*total_sum)
		return false;

	return true;
	
}


void SplicingGraph::refine_forward_trunk(ReadHash& read_hash, SuffixTree& right_tree) {

	if (node_set[0].sequence.length() <= g_read_length)
		return;
	const string& trunk = node_set[0].sequence;
	int trunk_len = trunk.length();
	int check_len = g_min_same_len*2;
	int read_id = -1;

	if (check_len + g_read_length > trunk_len)
		check_len = trunk_len - g_read_length;

	for (int i = 1; i <= check_len; ++i) {
		const string& read = trunk.substr(trunk_len-i-g_read_length, g_read_length);
		read_int_type read_int = get_read_int_type(read);
		read_id = read_hash.find_read(read_int);
		if (read_id  < 0)
			continue;
		const string& forward_trunk = forward_extend(right_tree, read_id);
		if (forward_trunk.length() >= g_read_length + g_min_exon_length) {
			node_set[0].sequence = trunk.substr(0, trunk_len-i-g_read_length) + forward_trunk;
			if (trunk_len-i-2*g_read_length+1 <= 0)
				set_reads_tag(read_hash, trunk, -4);
			else {
				const string& recover_str = trunk.substr(trunk_len-i-2*g_read_length+1);
				set_reads_tag(read_hash, recover_str, -4);
			}
			break;
		} else {
			set_reads_tag(read_hash,forward_trunk,-4);
		}
	}

}


void SplicingGraph::refine_reverse_trunk(ReadHash& read_hash, SuffixTree& left_tree) {

	if (node_set[0].sequence.length() <= g_read_length)
		return;
	const string& trunk = node_set[0].sequence;
	int trunk_len = trunk.length();
	int check_len = g_min_same_len*2;
	if (check_len + g_read_length > trunk_len)
		check_len = trunk_len - g_read_length;
	int read_id = 0;
	for (int i = 1; i <= check_len; ++i) {

		const string& read = trunk.substr(i, g_read_length);
		read_int_type read_int = get_read_int_type(read);
		read_id = read_hash.find_read(read_int);
		if (read_id < 0)
			continue;
		const string& reverse_trunk = reverse_extend(left_tree, read_id);
		if (reverse_trunk.length() >= g_read_length + g_min_exon_length) {
			node_set[0].sequence = reverse_trunk + trunk.substr(i+g_read_length);
			if (i + 2*g_read_length-1 > trunk_len)
				set_reads_tag(read_hash, trunk, -4);
			else {
				const string& recover_str = trunk.substr(0, i+g_read_length);
				set_reads_tag(read_hash, recover_str, -4);
			}
			break;
		} else {
			set_reads_tag(read_hash,reverse_trunk,-4);
		}
	}

}


void SplicingGraph::set_reads_pos_in_node(ReadHash& read_hash, int p) {

	if (p < 0 || p >= node_sum )
		return;
	if (node_set[p].sequence.length() < g_read_length)
		return;

	const string& sequence = node_set[p].sequence;
/*
	if (node_set[p].cov_reads.size() > 0) {
		for (int i = 0; i < node_set[p].cov_reads.size(); ++i)
			data_tag[node_set[p].cov_reads[i].first] = -4;
		node_set[p].cov_reads.clear();
	}
*/
	node_set[p].cov_reads.clear();
	read_int_type read_int;
	int read_id;
	pair<int, int> r_p;
	for (int i = 0; i <= sequence.length() - g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_id = read_hash.find_read(read_int);
		if (read_id < 0)
			continue;
		data_tag[read_id] = i;
		r_p.first = read_id;
		r_p.second = i;
		node_set[p].cov_reads.push_back(r_p);		
	}

}


bool SplicingGraph::refine_reverse_by_pair(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree, int p) {

//cout << "refine reverse by pair ..." << endl;
	if (node_set[p].sequence.length() < g_read_length+g_refine_check_len)
		return false;

	bool return_tag = false;
	vector<int> check_reads;
	int extend_time = 0;
	bool is_add = false;
	bool add_tag = true;
	int read_id, mate_id;

	while (1) {		
		const string& check_seq = node_set[p].sequence.substr(0, g_read_length+g_refine_check_len);
		check_reads.clear();
		set_reads_tag(read_hash, check_seq,check_reads, p);
		if (check_reads.size() == 0)
			break;
		is_add = false;
		for (int i = 0; i < check_reads.size(); ++i){
			add_tag = true;
			read_id = check_reads[i];
			if (read_id < g_divide_right_pos)
				continue;
			for (int j = 0; j < data_pair[read_id].size(); ++j) {
				if (data_tag[data_pair[read_id][j]] > -2)
					add_tag = false;
			}
			if (add_tag) {
				++extend_time;
				if (extend_time > 20)
					break;
				for (int j = 0; j < data_pair[read_id].size(); ++j) {
					mate_id = data_pair[read_id][j];
					string extend_seq = data[mate_id];
					data_tag[mate_id] = -1;
					const string& stop_seq = node_set[p].sequence.substr(0, g_min_same_len-1);
					is_add = forward_extend(right_tree, extend_seq, stop_seq, mate_id);
					if (is_add) {
						return_tag = true;
						const string& seq = reverse_extend(left_tree, mate_id);
						if (seq.length() > g_read_length)
							extend_seq = seq.substr(0,seq.length()-g_read_length) + extend_seq;
						node_set[p].sequence = extend_seq + node_set[p].sequence;
						break;
					} else {
						if (extend_seq.length() > 1000) {
							Node node1;
							node1.sequence = extend_seq;
							int n1 = add_node(node1);
							set_reads_tag(read_hash, extend_seq, n1);
						} else {
							set_reads_tag(read_hash,extend_seq,-4);
						}
					}
				}
			}
			if (is_add)
				break;
		}
		if (!is_add)
			break;
	}

	return return_tag;

}



bool SplicingGraph::refine_forward_by_pair(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree, int p) {

//cout << "refine forward by pair..." << endl;
	if (node_set[p].sequence.length() < g_read_length+g_refine_check_len)
		return false;

	bool return_tag = false;
	vector<int> check_reads;
	int extend_time = 0;
	bool is_add = false;
	bool add_tag = true;
	int read_id, mate_id;

	while (1) {		
		const string& check_seq = node_set[p].sequence.substr(node_set[p].sequence.length()-g_read_length-g_refine_check_len);
		check_reads.clear();
		set_reads_tag(read_hash, check_seq,check_reads, p);
		if (check_reads.size() == 0)
			break;
		bool is_add = false;
		for (int i = 0; i < check_reads.size(); ++i){
			add_tag = true;
			read_id = check_reads[i];
			if (read_id >= g_divide_right_pos)
				continue;
			for (int j = 0; j < data_pair[read_id].size(); ++j) {
				if (data_tag[data_pair[read_id][j]] > -2)
					add_tag = false;
			}
			if (add_tag) {
				++extend_time;
				if (extend_time > 20)
					break;
				for (int j = 0; j < data_pair[read_id].size(); ++j) {
					mate_id = data_pair[read_id][j];
					string extend_seq = data[mate_id];
					data_tag[mate_id] = -1;
					const string& stop_seq = node_set[p].sequence.substr(node_set[p].sequence.length()-g_min_same_len+1);
					is_add = reverse_extend(left_tree, extend_seq, stop_seq, mate_id);
					if (is_add) {
						return_tag = true;
						const string& seq = forward_extend(right_tree, mate_id);
						if (seq.length() > g_read_length)
							extend_seq = extend_seq + seq.substr(g_read_length);
						node_set[p].sequence = node_set[p].sequence + extend_seq;
						break;
					} else {
						if (extend_seq.length() > 1000) {
							Node node1;
							node1.sequence = extend_seq;
							int n1 = add_node(node1);
							set_reads_tag(read_hash, extend_seq, n1);
						} else {
							set_reads_tag(read_hash,extend_seq,-4);
						}
					}
				}
			}
			if (is_add)
				break;
		}
		if (!is_add)
			break;
	}

	return return_tag;

}


bool SplicingGraph::forward_extend(SuffixTree& right_tree, string& contig, const string& stop_seq, int seed) {

//cout << "forward extend..." << endl;
	if (contig.length() < g_read_length)
		return false;

	kmer_int_type kmer_int;
	trie_node* index;
	vector<vector<int> > candi_set;
	vector<int> candi_len;
	int overlap_len, pos;
	vector<int> read_vec;
	int read_id, mate_id, candi_id;
	bool is_extend = true;
	int contig_cov = data_cov[seed];
	int candi_cov = contig_cov*10000;

	while (is_extend) {

		const string& kmer = contig.substr(contig.length()-g_min_same_len);
		kmer_int = kmer_to_int(kmer);
		index = right_tree[kmer_int];
		if (index == NULL)
			break;
		candi_set.clear();
		candi_len.clear();
		overlap_len = g_min_same_len;
		while (index != NULL) {
			if ( (index->read_set).size() > 0){
				read_vec = index -> read_set;
				candi_set.push_back(read_vec);
				candi_len.push_back(overlap_len);
			}
			++overlap_len;
			pos = base_to_int(contig[contig.length()-overlap_len]);
			index = index -> children[pos];
		}		
	
		if (candi_len.size() == 0)
			break;

		is_extend = false;

		if (g_is_paired_end) {
			for (int i = candi_set.size()-1; i >= 0; --i) {
				read_vec = candi_set[i];
				for (int j = 0; j < read_vec.size(); ++j) {
					read_id = read_vec[j];
					int compatible_len = get_compatible_len(stop_seq, data[read_id].substr(g_read_length-g_min_same_len+1));
					if (compatible_len != 0){
						contig = contig.substr(0, contig.length()-candi_len[i]) + data[read_id].substr(0, g_read_length-compatible_len);
						data_tag[read_id] = -1;
						return true;
					}
					if (data_tag[read_id] > -2)
						continue;
					for (int k = 0; k < data_pair[read_id].size(); ++k) {
						mate_id = data_pair[read_id][k];
						if (data_tag[mate_id] > -2) {
							contig = contig + data[read_id].substr(candi_len[i]);
							data_tag[read_id] = -1;
							contig_cov = data_cov[read_id];
							is_extend = true;
							break;
						}
					}
					if (is_extend)
						break;
				}
				if(is_extend)
					break;				
			}

		}//if (g_is_paired_end)
		
		if (!is_extend) {
			for (int i = candi_set.size()-1; i >= 0; --i) {
				read_vec = candi_set[i];
				candi_id = -1;
				for (int j = 0; j < read_vec.size(); ++j) {
					read_id = read_vec[j];
					if (data_tag[read_id] > -2)
						continue;
					if (j == 0 || candi_cov < abs(data_cov[read_id] - contig_cov) ) {
						candi_id = read_id;
						candi_cov = abs(data_cov[read_id] - contig_cov);
					}
				}
				if (candi_id >= 0){
					contig = contig + data[candi_id].substr(candi_len[i]);
					is_extend = true;
					data_tag[candi_id] = -1;
					contig_cov = data_cov[candi_id];
					break;
				}				
			}
			
		} //if (!is_extend)		
	}

	return false;
}





bool SplicingGraph::reverse_extend(SuffixTree& left_tree, string& contig, const string& stop_seq, int seed) {

//cout << "reverse extend..." << endl;
	if (contig.length() < g_read_length)
		return false;

	kmer_int_type kmer_int;
	trie_node* index;
	vector<vector<int> > candi_set;
	vector<int> candi_len;
	int overlap_len, pos;
	vector<int> read_vec;
	int read_id, mate_id,candi_id;
	bool is_extend = true;
	int contig_cov = data_cov[seed];
	int candi_cov = contig_cov*10000;

	while (is_extend) {

		const string& kmer = contig.substr(0, g_min_same_len);
		kmer_int = kmer_to_int(kmer);
		index = left_tree[kmer_int];
		if (index == NULL)
			break;
		candi_set.clear();
		candi_len.clear();
		overlap_len = g_min_same_len;
		while (index != NULL) {
			if ( (index->read_set).size() > 0){
				read_vec = index -> read_set;
				candi_set.push_back(read_vec);
				candi_len.push_back(overlap_len);
			}
			pos = base_to_int(contig[overlap_len]);
			index = index -> children[pos];
			++overlap_len;
		}	
	
		if (candi_len.size() == 0)
			break;

		is_extend = false;

		if (g_is_paired_end) {
			for (int i = candi_set.size()-1; i >= 0; --i) {
				read_vec = candi_set[i];
				for (int j = 0; j < read_vec.size(); ++j) {
					read_id = read_vec[j];
					int compatible_len = get_compatible_len(stop_seq, data[read_id].substr(0,g_min_same_len-1));
					if (compatible_len != 0){
						contig = data[read_id].substr(compatible_len) + contig.substr(candi_len[i]);
						data_tag[read_id] = -1;
						return true;
					}
					if (data_tag[read_id] > -2)
						continue;
					for (int k = 0; k < data_pair[read_id].size(); ++k) {
						mate_id = data_pair[read_id][k];
						if (data_tag[mate_id] > -2) {
							contig = data[read_id].substr(0, g_read_length-candi_len[i]) + contig;
							data_tag[read_id] = -1;
							contig_cov = data_cov[read_id];
							is_extend = true;
							break;
						}
					}
					if (is_extend)
						break;
				}
				if(is_extend)
					break;				
			}

		}//if (g_is_paired_end)
		
		if (!is_extend) {
			for (int i = candi_set.size()-1; i >= 0; --i) {
				read_vec = candi_set[i];
				candi_id = -1;
				for (int j = 0; j < read_vec.size(); ++j) {
					read_id = read_vec[j];
					if (data_tag[read_id] > -2)
						continue;
					if (j == 0 || candi_cov < abs(data_cov[read_id] - contig_cov) ) {
						candi_id = read_id;
						candi_cov = abs(data_cov[read_id] - contig_cov);
					}
				}
				if (candi_id >= 0){
					contig = data[candi_id].substr(0, g_read_length-candi_len[i]) + contig;
					is_extend = true;
					data_tag[candi_id] = -1;
					contig_cov = data_tag[candi_id];
					break;
				}				
			}
		} //if (!is_extend)		
	}
	
	return false;
}


void SplicingGraph::branch_extend_by_coverage(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree) {

	int current_node_sum = node_sum;
	for (int i = 0; i < current_node_sum; ++i)
		forward_extend_by_coverage(read_hash, right_tree, i);
	current_node_sum = node_sum;
	for (int i = 0; i < current_node_sum; ++i)
		reverse_extend_by_coverage(read_hash, left_tree, i);

}


void SplicingGraph::forward_extend_by_coverage(ReadHash& read_hash, SuffixTree& right_tree, int p) {

	if (p < 0 || p >= node_sum)
		return;

	//cout << "check forward branch..." << endl;

	set_forward_check_pos(read_hash, p);

	pair<int, int> right_id;
	read_int_type read_int;
	int read_id;

	while (node_set[p].forward_check_pos.size() > 0) {

		pair<int, int> check_pos = node_set[p].forward_check_pos.front();
		node_set[p].forward_check_pos.pop_front();
		if (node_set[p].sequence.length() < check_pos.second + g_read_length)
			break;
		if (check_pos.first <= g_read_length)
			break;

		for (int i = check_pos.first; i < check_pos.second; ++i) {
			const string& seed = node_set[p].sequence.substr(i-g_read_length, g_read_length);
			read_int = get_read_int_type(seed);
			read_id = read_hash.find_read(read_int);
			if (read_id < 0)
				continue;
			right_id.first = -1;
			right_id.second = -1;
			const string& forward_seq = forward_extend(right_tree, read_id);
			set_reads_tag(read_hash, forward_seq, -4);

			if (forward_seq.length() <= g_read_length + 10  || (g_is_paired_end &&(!pair_support(read_hash, forward_seq,p)))) 
				continue;
			const string& reach_seq = forward_seq.substr(forward_seq.length()-g_min_same_len);
			string::size_type start2 = node_set[p].sequence.substr(i).find(reach_seq);
			int right_node = -1;
			if (start2 == string::npos) {
				for (int k = 0; k < node_sum; ++k) {
					if (k ==p)
						continue;
					start2 = node_set[k].sequence.find(reach_seq);
					if (start2 != string::npos) {
						if (g_is_paired_end && (!pair_support(read_hash, forward_seq,k)))
							continue;
						right_id.first = k;
						right_id.second = start2;
						right_node = k;
						break;
					}
				}

			} else {
				right_id.first = p;
				right_node = p;
				right_id.second = start2;
			}


			if (right_id.second == -1 && forward_seq.length() < g_min_exon_length + g_read_length)
				continue;

			if (right_id.second == -1) {//add branch
				Node node2 = node_set[p].sequence.substr(0, i);
				int n2 = add_node(node2);
				node_set[p].sequence = node_set[p].sequence.substr(i);
				for (int j = 0; j < node_set[p].parents.size(); ++j) {
					int p_i = node_set[p].parents[j];
					node_set[p_i].delete_child(p);
					node_set[p_i].add_child(n2);
				}
				node_set[n2].parents = node_set[p].parents;
				node_set[p].parents.clear();
				node_set[p].add_parent(n2);
				node_set[n2].add_child(p);
				if (node_set[p].sequence.length() < g_min_exon_length && node_set[p].children.size() == 0) {
					set_reads_tag(read_hash, node_set[p].sequence, -4);
					node_set[p].sequence = forward_seq.substr(g_read_length);
					set_reads_tag(read_hash, forward_seq, p);
					forward_extend_by_coverage(read_hash, right_tree, p);
					return;
				} else {
					Node node1;
					node1.sequence = forward_seq.substr(g_read_length);
					int n1 = add_node(node1);
					node_set[n2].add_child(n1);
					node_set[n1].add_parent(n2);
					set_reads_tag(read_hash, forward_seq, n1);
					forward_extend_by_coverage(read_hash, right_tree, p);
					return;
				}
	
			} else { //add bubble
				bool fail_add = false;
				int length = forward_seq.length() - g_read_length - right_id.second;
				int start1 = i-g_read_length;
				string bubble_seq;
				if (length > 0)
					bubble_seq = forward_seq.substr(g_read_length, length);
				if ( p == right_node) {
					if ((!fail_add) && start2 <= start1)
						fail_add = true;
					int distance = 0;
					if (!fail_add)
						distance = start2 - start1 - g_read_length;
					if (distance <= 0)
						fail_add = true;
					if (distance + length <= 4)
						fail_add = true;
					if (!fail_add) {
						Node node1, node2, node3;
						int n1 = -1;
						int n2 = -1;
						int n3 = -1;
						if (length > 0) {
							node1.sequence = bubble_seq;
							n1 = add_node(node1);
							set_reads_tag(read_hash, forward_seq, n1);
						}
						node2.sequence = node_set[p].sequence.substr(0, start1+g_read_length);
						node3.sequence = node_set[p].sequence.substr(start1+g_read_length, distance);
						node_set[p].sequence = node_set[p].sequence.substr(start2);
						n2 = add_node(node2);
						n3 = add_node(node3);
						for (int j = 0; j < node_set[p].parents.size(); ++j) {
							int p_i = node_set[p].parents[j];
							node_set[p_i].delete_child(p);
							node_set[p_i].add_child(n2);
						}
						node_set[n2].parents = node_set[p].parents;
						node_set[p].parents.clear();
						node_set[n2].add_child(n3);
						node_set[n3].add_parent(n2);
						node_set[n3].add_child(p);
						node_set[p].add_parent(n3);
						if (length > 0) {
							node_set[n1].add_child(p);
							node_set[p].add_parent(n1);
							node_set[n2].add_child(n1);
							node_set[n1].add_parent(n2);
						} else {
							node_set[n2].add_child(p);
							node_set[p].add_parent(n2);
						}
						forward_extend_by_coverage(read_hash, right_tree, n3);
						forward_extend_by_coverage(read_hash, right_tree, p);
						return;
					}

				} else { //p != right_node
					set<int> checked;
					if (right_node == 0 || is_circle(p, right_node, checked))
						fail_add = true;							
					if (!fail_add) {
						Node node1, node2, node3;
						int n1 = -1;
						int n2 = -1;
						int n3 = -1;
						if (length > 0) {
							node1.sequence = bubble_seq;
							n1 = add_node(node1);
							set_reads_tag(read_hash, forward_seq, n1);
						}
						if (start1+g_read_length < node_set[p].sequence.length()) {
							node2.sequence = node_set[p].sequence.substr(0, start1+g_read_length);
							n2 = add_node(node2);
							node_set[p].sequence = node_set[p].sequence.substr(start1+g_read_length);
							for (int j = 0; j < node_set[p].parents.size(); ++j) {
								int p_i = node_set[p].parents[j];
								node_set[p_i].delete_child(p);
								node_set[p_i].add_child(n2);
							}
							node_set[n2].parents = node_set[p].parents;
							node_set[p].parents.clear();
							node_set[n2].add_child(p);
							node_set[p].add_parent(n2);
							if (length > 0) {
								node_set[n2].add_child(n1);
								node_set[n1].add_parent(n2);
							}
						} else {
							if (length > 0) {
								node_set[p].add_child(n1);
								node_set[n1].add_parent(p);
							}
						}
						if (start2 > 0) {
							node3.sequence = node_set[right_node].sequence.substr(0, start2);
							n3 = add_node(node3);
							node_set[right_node].sequence = node_set[right_node].sequence.substr(start2);
							for (int j = 0; j < node_set[right_node].parents.size(); ++j) {
								int p_i = node_set[right_node].parents[j];
								node_set[p_i].delete_child(right_node);
								node_set[p_i].add_child(n3);
							}
							node_set[n3].parents = node_set[right_node].parents;
							node_set[right_node].parents.clear();
							node_set[n3].add_child(right_node);
							node_set[right_node].add_parent(n3);									
						}
						if (length > 0) {
							node_set[right_node].add_parent(n1);
							node_set[n1].add_child(right_node);
						} else {
							if (start1+g_read_length < node_set[p].sequence.length()) {
								node_set[n2].add_child(right_node);
								node_set[right_node].add_parent(n2);
							} else {
								node_set[p].add_child(right_node);
								node_set[right_node].add_parent(p);
							}
						}
						forward_extend_by_coverage(read_hash, right_tree, p);
						return;
					}
				}


				if (fail_add && forward_seq.length() > g_min_exon_length + g_read_length){
						
					Node node2 = node_set[p].sequence.substr(0, i);
					int n2 = add_node(node2);
					node_set[p].sequence = node_set[p].sequence.substr(i);
					for (int j = 0; j < node_set[p].parents.size(); ++j) {
						int p_i = node_set[p].parents[j];
						node_set[p_i].delete_child(p);
						node_set[p_i].add_child(n2);
					}
					node_set[n2].parents = node_set[p].parents;
					node_set[p].parents.clear();
					node_set[p].add_parent(n2);
					node_set[n2].add_child(p);
					if (node_set[p].sequence.length() < g_min_exon_length && node_set[p].children.size() == 0) {
						set_reads_tag(read_hash, node_set[p].sequence, -4);
						node_set[p].sequence = forward_seq.substr(g_read_length);
						forward_extend_by_coverage(read_hash, right_tree, p);
						return;
					} else {
						Node node1;
						node1.sequence = forward_seq.substr(g_read_length);
						int n1 = add_node(node1);
						set_reads_tag(read_hash, forward_seq, n1);
						node_set[n2].add_child(n1);
						node_set[n1].add_parent(n2);							
						forward_extend_by_coverage(read_hash, right_tree, p);
						return;
					}

				} //if (fail_add)
			}//add bubble

		} // for

	}

}



void SplicingGraph::reverse_extend_by_coverage(ReadHash& read_hash, SuffixTree& left_tree, int p) {

	if (p < 0 || p >= node_sum)
		return;

	set_reverse_check_pos(read_hash, p);
//cout << "reverse  check " << p << endl;

	pair<int, int> right_id;
	read_int_type read_int;
	int read_id;

	while (node_set[p].reverse_check_pos.size() > 0) {

		pair<int, int> check_pos = node_set[p].reverse_check_pos.front();
		node_set[p].reverse_check_pos.pop_front();
		if (node_set[p].sequence.length() <= check_pos.second + g_read_length)
			break;
		if (check_pos.first <= 0)
			break;
//cout << check_pos.first << " " << check_pos.second <<" "<< node_set[p].sequence.length() << endl;
		for (int i = check_pos.first; i < check_pos.second; ++i) {
			const string& seed = node_set[p].sequence.substr(i, g_read_length);
			read_int = get_read_int_type(seed);
			read_id = read_hash.find_read(read_int);
			if (read_id < 0)
				continue;
			right_id.first = -1;
			right_id.second = -1;
			const string& reverse_seq = reverse_extend(left_tree, read_id);
			set_reads_tag(read_hash, reverse_seq, -4);

			string::size_type start2;
			int right_node = -1;
			const string& reach_seq = reverse_seq.substr(0, g_min_same_len);
			for (int k = 0; k < node_sum; ++k) {
				if (k ==p)
					continue;
				start2 = node_set[k].sequence.find(reach_seq);
				if (start2 != string::npos) {
					if (g_is_paired_end && (!pair_support(read_hash, reverse_seq,k)))
						continue;
					right_id.first = k;
					right_id.second = start2;
					right_node = k;
					break;
				}
			}

			if ((g_is_paired_end && (!pair_support(read_hash, reverse_seq, p))) || reverse_seq.length() <= g_read_length + 10) 
				continue;
			
			if (right_id.second == -1 && reverse_seq.length() < g_min_exon_length + g_read_length)
				continue;

			if (right_id.second == -1) {//add branch
				Node node2 = node_set[p].sequence.substr(0, i);
				int n2 = add_node(node2);
				node_set[p].sequence = node_set[p].sequence.substr(i);
				for (int j = 0; j < node_set[p].parents.size(); ++j) {
					int p_i = node_set[p].parents[j];
					node_set[p_i].delete_child(p);
					node_set[p_i].add_child(n2);
				}
				node_set[n2].parents = node_set[p].parents;
				node_set[p].parents.clear();
				node_set[p].add_parent(n2);
				node_set[n2].add_child(p);
				if (node_set[n2].sequence.length() < g_min_exon_length && node_set[n2].parents.size() == 0) {
					set_reads_tag(read_hash, node_set[n2].sequence, -4);
					node_set[n2].sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
					set_reads_tag(read_hash, reverse_seq, n2);
					reverse_extend_by_coverage(read_hash, left_tree, p);
					return;
				} else {
					Node node1;
					node1.sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
					int n1 = add_node(node1);
					node_set[n1].add_child(p);
					node_set[p].add_parent(n1);
					set_reads_tag(read_hash, reverse_seq, n1);
					reverse_extend_by_coverage(read_hash, left_tree, p);
					return;
				}
	
			} else { //add bubble
				bool fail_add = false;
				int length = reverse_seq.length() - g_read_length - right_id.second;
				int start1 = i;
				string bubble_seq;
				if (length > 0)
					bubble_seq = reverse_seq.substr(right_id.second, length);
				set<int> checked;
				if (right_node == 0 || is_circle(right_node, p, checked))
					fail_add = true;							
				if (!fail_add) {
					Node node1, node2, node3;
					int n1 = -1;
					int n2 = -1;
					int n3 = -1;
					if (length > 0) {
						node1.sequence = bubble_seq;
						n1 = add_node(node1);
						set_reads_tag(read_hash, reverse_seq, n1);
					}
					if (start1 > 0) {
						node2.sequence = node_set[p].sequence.substr(0, start1);
						n2 = add_node(node2);
						node_set[p].sequence = node_set[p].sequence.substr(start1);
						for (int j = 0; j < node_set[p].parents.size(); ++j) {
							int p_i = node_set[p].parents[j];
							node_set[p_i].delete_child(p);
							node_set[p_i].add_child(n2);
						}
						node_set[n2].parents = node_set[p].parents;
						node_set[p].parents.clear();
						node_set[n2].add_child(p);
						node_set[p].add_parent(n2);
						
					} 
					if (length > 0) {
						node_set[n1].add_child(p);
						node_set[p].add_parent(n1);
					}
					
					if (start2 + right_id.second < node_set[right_node].sequence.length()) {
						node3.sequence = node_set[right_node].sequence.substr(0, start2+right_id.second);
						n3 = add_node(node3);
						node_set[right_node].sequence = node_set[right_node].sequence.substr(start2 + right_id.second);
						for (int j = 0; j < node_set[right_node].parents.size(); ++j) {
							int p_i = node_set[right_node].parents[j];
							node_set[p_i].delete_child(right_node);
							node_set[p_i].add_child(n3);
						}
						node_set[n3].parents = node_set[right_node].parents;
						node_set[right_node].parents.clear();
						node_set[n3].add_child(right_node);
						node_set[right_node].add_parent(n3);									
					}
					if (length > 0) {
						if (start2 + right_id.second < node_set[right_node].sequence.length()) {
							node_set[n1].add_parent(n3);
							node_set[n3].add_child(n1);
						} else {
							node_set[n1].add_parent(right_node);
							node_set[right_node].add_child(n1);
						}
					} else {
						if (start2 + right_id.second < node_set[right_node].sequence.length()) {
							node_set[n3].add_child(p);
							node_set[p].add_parent(n3);
						} else {
							node_set[right_node].add_child(p);
							node_set[p].add_parent(right_node);
						}
					}
					reverse_extend_by_coverage(read_hash, left_tree, p);
					return;
				}

				if (fail_add && reverse_seq.length() > g_min_exon_length + g_read_length){
					Node node1;
					int n1 = -1;
					if (i == 0) {
						if (node_set[p].parents.size() == 0) {
							node_set[p].sequence = reverse_seq + node_set[p].sequence.substr(g_read_length);
							set_reads_tag(read_hash, reverse_seq, p);
						} else {	
							node1.sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
							n1 = add_node(node1);
							set_reads_tag(read_hash, reverse_seq, n1);
							node_set[n1].add_child(p);
							node_set[p].add_parent(n1);							
						}
					} else {
						Node node2 = node_set[p].sequence.substr(0, i);
						int n2 = add_node(node2);
						node_set[p].sequence = node_set[p].sequence.substr(i);
						for (int j = 0; j < node_set[p].parents.size(); ++j) {
							int p_i = node_set[p].parents[j];
							node_set[p_i].delete_child(p);
							node_set[p_i].add_child(n2);
						}
						node_set[n2].parents = node_set[p].parents;
						node_set[p].parents.clear();
						node_set[p].add_parent(n2);
						node_set[n2].add_child(p);
						node1.sequence = reverse_seq.substr(0, reverse_seq.length()-g_read_length);
						n1 = add_node(node1);
						set_reads_tag(read_hash, reverse_seq, n1);
						node_set[n1].add_child(p);
						node_set[p].add_parent(n1);
					}							
					reverse_extend_by_coverage(read_hash, left_tree, p);
					return;

				} //if (fail_add)
			}//add bubble

		} // for

	}

}


void SplicingGraph::set_forward_check_pos(ReadHash& read_hash, int p) {

//cout << "set check" << endl;
	node_set[p].forward_check_pos.clear();
	if (node_set[p].sequence.length() <= g_read_length + 1)
		return;

	set_reads_pos_in_node(read_hash, p);

	if (node_set[p].cov_reads.size() < 10)
		return;

	int node_len = node_set[p].sequence.length();
	vector<float> node_cov(node_len, 0.0);
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i)
		for (int j = node_set[p].cov_reads[i].second; j < node_set[p].cov_reads[i].second+g_read_length; ++j)
			node_cov[j] = node_cov[j] + 1.0;
///*
	int start, len, end;
	pair<int, int> fp;
	len = 0;
	start = 0;
	end = 0;

	for (int i = g_read_length; i < (int)node_cov.size()-1; ++i) {

		if (node_cov[i] <= node_cov[i+1]){

			if (len >= 10 && node_cov[end] < 0.5*node_cov[start]) {

				fp.first = start-g_read_length;
				fp.second = start-g_read_length+10;
				node_set[p].forward_check_pos.push_back(fp);
			}
			start = 0;
			end = 0;
			len = 0;
		
		} else {

			if (start > 0) {
				++len;
				end = i+1;
			} else {
				start = i+1;
			}
		}
	}

	int id = -1;
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i) {
		id = node_set[p].cov_reads[i].first;
		data_tag[id] = p;
	}

}


void SplicingGraph::set_reverse_check_pos(ReadHash& read_hash, int p) {

//cout << "set check" << endl;
	node_set[p].reverse_check_pos.clear();

	if (node_set[p].sequence.length() <= g_read_length + 1)
		return;

	set_reads_pos_in_node(read_hash, p);

	if (node_set[p].cov_reads.size() < 10)
		return;

	int node_len = node_set[p].sequence.length();
	vector<float> node_cov(node_len, 0.0);
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i)
		for (int j = node_set[p].cov_reads[i].second; j < node_set[p].cov_reads[i].second+g_read_length; ++j)
			node_cov[j] = node_cov[j] + 1.0;
//*	
	int start, len, end;
	pair<int, int> fp;
	start = 0;
	end = 0;
	len = 0;
	for (int i = g_read_length; i < (int)node_cov.size()-1; ++i) {

		if (node_cov[i] >= node_cov[i+1]){
			if (len >= 10 && node_cov[end] > 0.5*node_cov[start]) {
				fp.first = end-10;
				fp.second = end;
				node_set[p].reverse_check_pos.push_back(fp);
			}
			start = 0;
			end = 0;
			len = 0;
		
		} else {

			if (start > 0) {
				++len;
				end = i+1;
			} else {
				start = i+1;
			}
		}
	}

	int id = -1;
	for (int i = 0; i < node_set[p].cov_reads.size(); ++i) {
		id = node_set[p].cov_reads[i].first;
		data_tag[id] = p;
	}

}




bool SplicingGraph::pair_support(ReadHash& read_hash, const string& sequence, int p) {

//cout << "pair support " << endl;
	if (p < 0 || p >= node_sum)
		return false;

	if (!g_is_paired_end)
		return true;

	int support = 0;
	int support1 = 0;
	read_int_type read_int;
	int read_id;

	for (int i = 0; i < sequence.length()-g_read_length; ++i) {
		const string& read_seq = sequence.substr(i, g_read_length);
		read_int = get_read_int_type(read_seq);
		read_id = read_hash.find_read(read_int);
		if (read_id < 0)
			continue;
		for (int j = 0; j < data_pair[read_id].size(); ++j) {
			if (data_tag[data_pair[read_id][j]] == p) {
				++support;
				if (support >= g_min_reads_support_branch)
					return true;
			} else {
				if (data_tag[data_pair[read_id][j]] >= 0) {
					++support1;
					if (support1 >= g_min_reads_support_branch*3)
						return true;
				}
			}
		}
		
	}
	
	return false;

}



bool SplicingGraph::is_circle(int p, int q, set<int>& checked) {

	bool flag = false;
	if (node_set[q].children.empty() || checked.find(q) != checked.end())
		return false;
	else
		checked.insert(q);

	 vector<int>::iterator it;
	 for (it = node_set[q].children.begin(); it != node_set[q].children.end(); it++) {

		 if ((*it) == p) {
			 flag = true;
			 break;
		 } else {
			 flag = is_circle(p, *it, checked);
			 if (flag)
				 break;
		 }

	 }

	 return flag;

}


void SplicingGraph::refine_graph(ReadHash& read_hash, SuffixTree& left_tree, SuffixTree& right_tree) {

	//cout << "refine graph..." << endl;
	int num = 0;
	for (int i = 0; i < node_sum; ++i) {
		bool is_forward = false;
		bool is_reverse = false;
		if (node_set[i].children.size() == 0) {
			num = 0;
			is_forward=refine_forward_by_pair(read_hash, left_tree, right_tree, i);
		}
		if (node_set[i].parents.size() == 0) {
			num = 0;
			is_reverse=refine_reverse_by_pair(read_hash, left_tree, right_tree, i);
		}
		if (is_forward || is_reverse)
			set_reads_tag(read_hash, node_set[i].sequence, i);

	}

}


void SplicingGraph::topological_sort() {
	
	Node s,t;
	s.sequence = "start";
	t.sequence = "end";
	int s_id = add_node(s);
	int t_id = add_node(t);
	vector<int> node_color;

	for (int i = 0; i < node_sum; ++i) {

		if (node_set[i].parents.size() == 0 && i != s_id && i != t_id){//

			node_set[i].add_parent(s_id);
			node_set[s_id].add_child(i);

		} 

		if (node_set[i].children.size() == 0 && i != s_id && i != t_id) {//

			node_set[i].add_child(t_id);
			node_set[t_id].add_parent(i);

		}

		node_color.push_back(0);

	}

	dfs_visit(s_id, node_color);

}


void SplicingGraph::dfs_visit(int i, vector<int>& node_color) {

	node_color[i] = 1;

	if (node_set[i].children.size() == 0) 
		node_order.push_back(i);
	else {

		for (int j = 0 ; j < node_set[i].children.size(); ++j) {

			if (node_color[node_set[i].children[j]] == 0)
				dfs_visit(node_set[i].children[j], node_color);

		}

		node_order.push_back(i);

	}

}


void SplicingGraph::set_coverage_of_nodes(ReadHash& read_hash) {

	for (int i = 0; i < node_sum; ++i) {

		if (node_set[i].sequence.length() < g_read_length) {
			node_set[i].node_coverage = 0;
			continue;
		}
		//set_reads_pos_in_node(read_hash, i);
		node_set[i].node_coverage = compute_coverage(read_hash, node_set[i].sequence);

	}

}


float SplicingGraph::compute_coverage(ReadHash& read_hash, const string& sequence) {


	if (sequence.length() < g_read_length)
		return 0.0;
	float total = 0.0;

	int read_id;
	for (int i = 0; i <= sequence.length()-g_read_length; ++i) {
		const string& read = sequence.substr(i, g_read_length);
		read_int_type read_int =  get_read_int_type(read);
		read_id = read_hash.find_read(read_int);
		if (read_id < 0 || data_tag[read_id] == -2)
			continue;
		total = total + data_pair[read_id].size();
	}

	total = total * g_read_length / sequence.length();
	return total;
}



void SplicingGraph::trim_graph(ReadHash& read_hash) {

	set_coverage_of_nodes(read_hash);
	bool has_trimed = false;
	map<pair_t, double> edge_coverages;
	vector<pair_t> edges;
	get_coverage_of_edges(read_hash, edge_coverages, edges);
	map<int, double> total_out_coverages;
	map<int, double> total_in_coverages;


	for (int i = 0; i < edges.size(); ++i) {

		int source = edges[i].first;
		if (total_out_coverages.find(source) == total_out_coverages.end())
			total_out_coverages[source] = edge_coverages[edges[i]];
		else
			total_out_coverages[source] += edge_coverages[edges[i]];

		int target = edges[i].second;
		if (total_in_coverages.find(target) == total_in_coverages.end())
			total_in_coverages[target] = edge_coverages[edges[i]];
		else
			total_in_coverages[target] += edge_coverages[edges[i]];

	}

	for (int i = 0; i < edges.size(); ++i) {

		int source = edges[i].first;
		int target = edges[i].second;

		//if ((node_set[source].parents.empty() && node_set[source].children.size() == 1) || 
		//	(node_set[target].children.empty() && node_set[target].parents.size() == 1)) {

				const string& check_edge = get_edge_sequence(source, target);
				if ((int)check_edge.length() < g_read_length)
					continue;
				double e_cov = edge_coverages[edges[i]];
				double flanking_node_cov = node_set[source].node_coverage > node_set[target].node_coverage ? node_set[source].node_coverage : node_set[target].node_coverage;
				
				if (e_cov < g_min_junction_coverage ||
					e_cov < g_min_ratio_welds * flanking_node_cov ||
					e_cov < g_min_ratio_branch * total_out_coverages[source] ||
					e_cov < g_min_ratio_branch * total_in_coverages[target] ||
					(total_in_coverages.find(source) != total_in_coverages.end() && e_cov < g_min_ratio_in_out * total_in_coverages[source]) ||
					(total_out_coverages.find(target) != total_out_coverages.end() && e_cov < g_min_ratio_in_out * total_out_coverages[target])
					) {

						has_trimed = true;
						set_reads_tag(read_hash, check_edge, -4);
						node_set[source].delete_child(target);
						node_set[target].delete_parent(source);
				}

		//}

	}
	set_parents();
}

void SplicingGraph::get_coverage_of_edges(ReadHash& read_hash, map<pair_t, double>& edge_coverages, vector<pair_t>& edges) {

	edges.clear();
	edge_coverages.clear();

	for (int i = 0; i < node_sum; ++i) {

		if (node_set[i].children.empty())
			continue;

		for (int j = 0; j < (int)node_set[i].children.size(); ++j)
			edges.push_back(pair_t(i, node_set[i].children[j]));

	}

	for (int i = 0; i < (int)edges.size(); ++i) {

		int source = edges[i].first;
		int target = edges[i].second;
		const string& edge = get_edge_sequence(source, target);
		edge_coverages[edges[i]] = compute_coverage(read_hash,edge);

	}

}


void SplicingGraph::get_coverage_of_edges(ReadHash& read_hash, map<pair_t, double>& edge_coverages) {

	vector<pair_t> edges;

	for (int i = 0; i < node_sum; ++i) {

		if (node_set[i].children.empty())
			continue;

		for (int j = 0; j < (int)node_set[i].children.size(); ++j)
			edges.push_back(pair_t(i, node_set[i].children[j]));

	}

	for (int i = 0; i < (int)edges.size(); ++i) {

		int source = edges[i].first;
		int target = edges[i].second;
		const string& edge = get_edge_sequence(source, target);
		edge_coverages[edges[i]] = compute_coverage(read_hash, edge);

	}

}




string SplicingGraph::get_edge_sequence(int s, int t) {

	if (s < 0 || s >= node_sum)
		return "";
	if (t < 0 || t >= node_sum)
		return "";

	int length = g_read_length - 1;
	int start = (int)node_set[s].sequence.length() > length ? static_cast<int>(node_set[s].sequence.length())-length : 0;
	string edge = node_set[s].sequence.substr(start);

	while ((int)edge.length() < length) {

		if (node_set[s].parents.size() >= 1) {

			int parent = node_set[s].parents[0];
			int remain = length - edge.length();
			int start_p = (int)node_set[parent].sequence.length() > remain ? node_set[parent].sequence.length()-remain : 0;
			edge = node_set[parent].sequence.substr(start_p) + edge;
			s = parent;

		} else {

			break;

		}

	}
	int pre_len = edge.length();
	if ((int)node_set[t].sequence.length() < length)
		edge = edge + node_set[t].sequence;
	else
		edge = edge + node_set[t].sequence.substr(0, length);

	while ( (int)edge.length() < pre_len+length) {
		
		if (node_set[t].children.size() >= 1) {

			int child = node_set[t].children[0];
			int remain =pre_len + length - edge.length();
			if (node_set[child].sequence.length() >= remain)
				edge = edge + node_set[child].sequence.substr(0,remain);
			else 
				edge = edge + node_set[child].sequence;
			t = child;

		} else {

			break;

		}

	}

	return edge;

}


void SplicingGraph::get_transcripts_by_all(ReadHash& read_hash, vector<pair<string,float> >& transcripts) {

	//cout << "get_transcripts..." << endl;

	recover_reads();
	
	for (int i = 0; i < node_set[node_order[node_order.size()-1]].children.size(); ++i) {
		int root_child = node_set[node_order[node_order.size()-1]].children[i];
		vector<int> path_node;
		path_node.push_back(root_child);
		get_all_path(root_child, path_node,transcripts);

	}
	check_transcripts(read_hash, transcripts);
}


void SplicingGraph::get_all_path(int node, vector<int>& path_node, vector<pair<string,float> >& transcripts) {
	for (int i = 0; i < node_set[node].children.size(); ++i) {
		set<int> checked;
		if (node_set[node].children[i] == node_order[0] || is_circle(node, node_set[node].children[i], checked)) {
			string transcript;
			for (int j = 0; j < path_node.size(); ++j) {
				if (j == 0){
					transcript = node_set[path_node[j]].sequence;
				} else {
					transcript = transcript + node_set[path_node[j]].sequence;
				}
			}
			pair<string, float> new_trans;
			new_trans.first = transcript;
			new_trans.second = 2.0;	
			transcripts.push_back(new_trans);
		} else {
			path_node.push_back(node_set[node].children[i]);
			get_all_path(node_set[node].children[i], path_node, transcripts);
			path_node.pop_back();
		}	
	}
}



void SplicingGraph::get_transcripts_by_cov_and_len(ReadHash& read_hash, vector<pair<string,float> >& transcripts) {

	//cout << "get_transcripts..." << endl;
	set_coverage_of_nodes(read_hash);
	float min_x = 1000.0;
	recover_reads();
	for (int i = 0; i < node_sum; ++i) {
		if (min_x > node_set[i].node_coverage && node_set[i].node_coverage > 0) 
			min_x = node_set[i].node_coverage;
	}
	int while_tag = 0;	
	map<pair_t, int> edge_tag;
	for (int i = 0; i < node_sum; i++) {
		for (int j = 0; j < node_set[i].children.size(); j++) {
			pair_t edge_t;
			edge_t.first = i;
			edge_t.second = node_set[i].children[j];
			edge_tag[edge_t] = 1;
			++while_tag;
		}
	}

	for (int i = 0; i < node_sum; i++) {	
		if (node_set[i].parents.size()==1 && node_set[i].parents[0] == node_order[node_sum-1] && node_set[i].children.size() == 1 && node_set[i].children[0] == node_order[0]) {

			while_tag = while_tag -2;
			pair<string, float> new_trans;
			new_trans.first = node_set[i].sequence;
			new_trans.second = node_set[i].node_coverage;
			transcripts.push_back(new_trans);
			node_set[i].delete_child(node_order[0]);
			node_set[node_order[0]].delete_parent(i);
			node_set[i].delete_parent(node_order[node_sum-1]);
			node_set[node_order[node_sum-1]].delete_child(i);
		}
	}

	map<pair_t, double> edge_coverages;
	get_coverage_of_edges(read_hash, edge_coverages);

	vector<float> min_cov;
	for (int i = 0; i < node_sum; i++)
		min_cov.push_back(node_set[i].node_coverage);
	int tag_sum = 0;
	tag_sum = tag_sum + node_set[node_order[node_order.size()-1]].children.size();
	tag_sum = tag_sum + node_set[node_order[0]].parents.size();
	map<pair_t, int>::iterator its;
	while (while_tag > tag_sum) {

		vector<int> node_parents;
		for (int i = 0; i < node_sum; i++)
			node_parents.push_back(-1);
		for (int i = node_sum-2; i > 0; i--) {

			int v = node_order[i];
			if (node_set[v].parents.size() == 1 && node_set[v].parents[0] == node_order[node_sum-1])
				continue;
			int max_cov = 0;
			for (int j = 0; j < node_set[v].parents.size(); j++) {
				int u = node_set[v].parents[j];
				float edge_cov = edge_coverages[pair_t(u,v)];
				if (edge_cov == 0.0)
					continue;
				if (edge_cov > min_cov[u])
					edge_cov = min_cov[u];
				if (max_cov < edge_cov) {
					max_cov = edge_cov;
					node_parents[v] = u;
				}
			}		
			min_cov[v] = max_cov;

		}
		int max_leaf = node_set[node_order[0]].parents[0];
		float max_leaf_cov = min_cov[node_set[node_order[0]].parents[0]];

		int max_len = 0;
		for (int i = 0; i < node_set[node_order[0]].parents.size(); i++) {
			int v = node_set[node_order[0]].parents[i];
			int current_len = 0;
			while (node_parents[v] >= 0){
				v = node_parents[v];
				current_len = current_len + 1;
			}
			if (current_len > max_len) {
				max_len = current_len;
				max_leaf = node_set[node_order[0]].parents[i];
				max_leaf_cov = min_cov[node_set[node_order[0]].parents[i]];
			}
		}
				
		string trans_seq = node_set[max_leaf].sequence;
		float trans_cov = max_leaf_cov;
		int check_x = 0;
		while (1) {

			if (node_parents[max_leaf] >= 0) {
				edge_coverages[pair_t(node_parents[max_leaf], max_leaf)] = edge_coverages[pair_t(node_parents[max_leaf], max_leaf)] - max_leaf_cov;
				min_cov[max_leaf] = min_cov[max_leaf] - max_leaf_cov;
				if (min_cov[max_leaf] <= 0)
					min_cov[max_leaf] = min_x;
				check_x = check_x + edge_tag[pair_t(node_parents[max_leaf], max_leaf)];
				trans_seq = node_set[node_parents[max_leaf]].sequence + trans_seq;
				edge_tag[pair_t(node_parents[max_leaf], max_leaf)] = 0;
				max_leaf = node_parents[max_leaf];
			} else
				break;

		}
		if (trans_cov == 0.0)
			break;
		if (check_x == 0)
			break;
		pair<string, float> new_trans;
		new_trans.first = trans_seq;
		new_trans.second = trans_cov;
		transcripts.push_back(new_trans);
		while_tag = 0;
		for (its = edge_tag.begin(); its != edge_tag.end(); its++)
			while_tag = while_tag + its -> second;
	}

	check_transcripts(read_hash, transcripts);
}


void SplicingGraph::recover_reads() {

	for (int i = 0; i < node_sum; ++i) {
		
		for (int it = 0; it < node_set[i].cov_reads.size(); it++) {
			data_tag[node_set[i].cov_reads[it].first] = -4;
		}
	}
}


void SplicingGraph::check_transcripts(ReadHash& read_hash, vector<pair<string,float> >& transcripts) {

	//cout << "check_transcripts.." << endl;

	vector<int> add_trans;
	vector<int> set_reads_vec;
	vector<int> recover_reads_vec;
	set<int> used_reads;
	float trans_cov;
	int mate_id, read_id;
///*
	string max_transcript = "N";
	int max_tran_id = -1;
	for (int i = 0; i < transcripts.size(); ++i) 
		if (transcripts[i].first.length() > max_transcript.length()) {
			max_transcript = transcripts[i].first;
			max_tran_id = i;
		}

//*/
	vector<pair<int, int> > reads_pos;
	//reads_pos.reserve(10000);
	for (int i = 0; i < transcripts.size(); ++i) {
		string transcript = transcripts[i].first;
		reads_pos.clear();
		map_reads_to_sequence(read_hash, transcript, reads_pos);
		bool is_add = true;
///*
		if (is_seq_similar(max_transcript, transcript, 'F',0.05) && i != max_tran_id)
			is_add = false;
		if (is_seq_similar(max_transcript, transcript,'R', 0.05) && i != max_tran_id)
			is_add = false;

		trans_cov = reads_pos.size()*g_read_length* 1.0 / transcript.length();
		//if ((transcript.length() < 500 && trans_cov < g_min_trans_cov*4) || (transcript.length() < 800 && trans_cov < g_min_trans_cov*2) || (transcript.length() < 1000 && trans_cov < g_min_trans_cov*1.5))
			//is_add = false;
		if ((!is_add) || trans_cov < g_min_trans_cov || transcript.length() < g_min_transcript_length) {

			for (int j = 0; j < reads_pos.size(); ++j)
				data_tag[reads_pos[j].first] = -4;
			continue;
		}

		if (g_is_paired_end) {
			
			float pair_end_sum = 0.0;		
			vector<pair<int, int> > maybe_used;
			for (int j = 0; j < reads_pos.size(); ++j) {

				int read_id = reads_pos[j].first;
				for (int k = 0; k < data_pair[read_id].size(); ++k) {
					mate_id = data_pair[read_id][k];
					if (data_tag[read_id] >= 0 && data_tag[mate_id] >= 0) {
						maybe_used.push_back(reads_pos[j]);
						pair_end_sum = pair_end_sum + 1.0;
					} else {
						data_tag[read_id] = -4;
					}
				}
				
			}

			trans_cov = pair_end_sum * g_read_length / transcript.length();
//cout << pair_end_sum  << " " << reads_pos.size() << endl;
			if ((transcript.length() < 500 && trans_cov < g_min_trans_cov*4) ||
				(transcript.length() < 700 && trans_cov < g_min_trans_cov*2) || 				(transcript.length() < 1000 && trans_cov < g_min_trans_cov*1.2))
				is_add = false;
			if (is_add && pair_end_sum > 10) {
				int total_cov_len = 0;
				int current_cov_len = g_read_length;
				int current_cov_start = maybe_used[0].second;
				bool add_tag = true;
				
				for (int j = 0; j < maybe_used.size()-1; ++j) {
					if (maybe_used[j].second + g_read_length >= maybe_used[j+1].second) {
						add_tag = true;
						current_cov_len = maybe_used[j+1].second + g_read_length - current_cov_start;
					} else {
						add_tag = false;
						total_cov_len = total_cov_len + current_cov_len;
						current_cov_len = g_read_length;
						current_cov_start = maybe_used[j+1].second;
					}
				}
				if (add_tag = true)
					total_cov_len = total_cov_len + current_cov_len;
//cout << "total_cov_len " << total_cov_len << endl;
//cout << transcript.length() << endl;
				//if (total_cov_len > 0.8*transcript.length()){
				if ((total_cov_len > 0.8*transcript.length() && transcript.length() >1200) || (total_cov_len>0.9*transcript.length() && transcript.length() > 500) || total_cov_len > 0.99*transcript.length() || total_cov_len > 1500){
					transcripts[i].second = pair_end_sum * g_read_length / transcript.length();
					for (int j = 0; j < maybe_used.size(); ++j) {
						read_id = maybe_used[j].first;
						data_cov[read_id] = data_cov[read_id] -1;
						if (data_cov[read_id] < 0)
							data_tag[read_id] = -2;
						else
							data_tag[read_id] = -3;
					}
				} else {
					is_add = false;
				}
			} else { 
				is_add = false;
			}					
		}
		
		if (is_add) {
			add_trans.push_back(i);			
		}
		for (int j = 0; j < reads_pos.size(); ++j)
			data_tag[reads_pos[j].first] = -4;

	}

	vector<pair<string,float> > temp_trans;
	for (int i = 0; i < add_trans.size(); ++i) {
		int tid = add_trans[i];
		temp_trans.push_back(transcripts[tid]);
	}
	transcripts.clear();	
	transcripts = temp_trans;

}


void SplicingGraph::map_reads_to_sequence(ReadHash& read_hash, const string& sequence, vector<pair<int, int> >& reads_pos) {

	if (sequence.length() < g_read_length)
		return;

	read_int_type read_int;
	int read_id;
	pair<int, int> r_p;

	for (int i = 0; i <= sequence.length() - g_read_length; ++i) {
		const string& real_read = sequence.substr(i,g_read_length);
		read_int = get_read_int_type(real_read);
		read_id = read_hash.find_read(read_int);
		if (read_id < 0 || data_tag[read_id] == -2)
			continue;
		r_p.first = read_id;
		r_p.second = i;
		data_tag[read_id] = i;
		reads_pos.push_back(r_p);		
	}	
}




