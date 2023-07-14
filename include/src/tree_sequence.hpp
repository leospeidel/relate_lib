#ifndef TREE_SEQUENCE_HPP
#define TREE_SEQUENCE_HPP

#include <iostream>
#include <err.h>
#include <vector>
#include <random>
#include <map>
#include <cstring>

#include "gzstream.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "mutations.hpp"
#include "tskit.h"

#define check_tsk_error(val) if (val < 0) {\
	errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));\
}

//Struct storing on how many branches mutation maps
struct PropagateStructLocal{

	int num_carriers = 0;
	int num_flipped_carriers = 0;
	int best_branch = -1;
	int best_flipped_branch = -1;

};

struct PropagateStructGlobal{

	int num_correct_carriers, num_correct_noncarriers;
	int num_incorrect_carriers, num_incorrect_noncarriers;
	int best_branch, best_flipped_branch;
	int min, flipped_min;

};

void
PropagateMutationExact(Node& node, std::vector<int>& branches, std::vector<int>& branches_flipped, Leaves& sequences_carrying_mutations, PropagateStructLocal& report){

	if(node.child_left != NULL){

		PropagateStructLocal report_c1, report_c2;

		PropagateMutationExact(*node.child_left, branches, branches_flipped, sequences_carrying_mutations, report_c1);
		PropagateMutationExact(*node.child_right, branches, branches_flipped, sequences_carrying_mutations, report_c2);

		report.num_carriers = report_c1.num_carriers + report_c2.num_carriers;
		report.num_flipped_carriers = report_c1.num_flipped_carriers + report_c2.num_flipped_carriers; 
		float num_leaves = report.num_carriers + report.num_flipped_carriers;

		if(report.num_flipped_carriers/num_leaves < 0.03 && report_c1.best_branch != -1 && report_c2.best_branch != -1){
			if(report_c1.num_carriers > 0 && report_c2.num_carriers > 0){
				report.best_branch             = node.label;
			}else if(report_c1.num_carriers > 0){
				report.best_branch             = report_c1.best_branch;
			}else if(report_c2.num_carriers > 0){
				report.best_branch             = report_c2.best_branch;
			}else{
				assert(false);
			}
		}else{
			if(report_c1.best_branch != -1){
				branches.push_back(report_c1.best_branch);
			}
			if(report_c2.best_branch != -1){
				branches.push_back(report_c2.best_branch);
			}
			report.best_branch = -1;
		}

		if(report.num_carriers/num_leaves < 0.03 && report_c1.best_flipped_branch != -1 && report_c2.best_flipped_branch != -1){
			if(report_c1.num_flipped_carriers > 0 && report_c2.num_flipped_carriers > 0){
				report.best_flipped_branch             = node.label;
			}else if(report_c1.num_flipped_carriers > 0){
				report.best_flipped_branch             = report_c1.best_flipped_branch;
			}else if(report_c2.num_flipped_carriers > 0){
				report.best_flipped_branch             = report_c2.best_flipped_branch;
			}else{
				assert(false);
			}
		}else{
			if(report_c1.best_flipped_branch != -1){
				branches_flipped.push_back(report_c1.best_flipped_branch);
			}
			if(report_c2.best_flipped_branch != -1){
				branches_flipped.push_back(report_c2.best_flipped_branch);
			}
			report.best_flipped_branch = -1;
		}

	}else{

		if(sequences_carrying_mutations.member[node.label] == 1){
			report.num_carriers         = 1;
			report.num_flipped_carriers = 0;
			report.best_branch          = node.label;
			report.best_flipped_branch  = -1;
		}else{
			report.num_carriers         = 0;
			report.num_flipped_carriers = 1;
			report.best_flipped_branch  = node.label;
			report.best_branch          = -1;
		}

	}

}

int 
MapMutationExact(Tree& tree, Leaves& sequences_carrying_mutations, Muts::iterator it_mut){

	//if(sequences_carrying_mutations.num_leaves == 0 || sequences_carrying_mutations.num_leaves == N) return 1;
	if(sequences_carrying_mutations.num_leaves == 0) return 1;

	//I want to place the mutation on all branches necessary for no loss of information
	//start with all leaves
	//propagate up and count number of nodes needed.
	//choose flipped or non-flipped depending on which is less.
	std::vector<int> branches;
	std::vector<int> branches_flipped; 

	PropagateStructLocal report;
	PropagateMutationExact(*std::prev(tree.nodes.end(), 1), branches, branches_flipped, sequences_carrying_mutations, report);

	if( branches_flipped.size() == 0 ){

		assert(branches.size() > 0);
		(*it_mut).branch  = branches;
		for(std::vector<int>::iterator it = branches.begin(); it != branches.end(); it++){
			tree.nodes[*it].num_events += 1.0;
		}
		return branches.size();

	}else{

		if( branches.size() <= branches_flipped.size() && branches.size() > 0 ){ 

			(*it_mut).branch  = branches;
			for(std::vector<int>::iterator it = branches.begin(); it != branches.end(); it++){
				tree.nodes[*it].num_events += 1.0;
			}
			return branches.size();

		}else{

			(*it_mut).flipped = true;
			(*it_mut).branch  = branches_flipped;
			for(std::vector<int>::iterator it = branches_flipped.begin(); it != branches_flipped.end(); it++){
				tree.nodes[*it].num_events += 1.0;
			}
			return branches_flipped.size();

		}

	}

}

void
PropagateMutation(Node& node, int N, Leaves& sequences_carrying_mutations, PropagateStructGlobal& report){

	float total_carriers    = sequences_carrying_mutations.num_leaves; 
	float total_noncarriers = N - total_carriers;  //slightly inefficient 

	if(node.child_left != NULL){

		PropagateStructGlobal report2;

		PropagateMutation(*node.child_left, N, sequences_carrying_mutations, report);
		PropagateMutation(*node.child_right, N, sequences_carrying_mutations, report2);

		report.num_correct_carriers      += report2.num_correct_carriers;
		report.num_incorrect_noncarriers += report2.num_incorrect_noncarriers;
		report.num_incorrect_carriers     = total_carriers - report.num_correct_carriers;
		report.num_correct_noncarriers    = total_noncarriers - report.num_incorrect_noncarriers;

		int sum = report.num_incorrect_carriers + report.num_incorrect_noncarriers;

		//std::cerr << report.num_incorrect_carriers/total_carriers << std::endl; 
		bool necessary_condition = (((float) report.num_incorrect_carriers)/total_carriers < 0.3);
		necessary_condition *= (((float) report.num_incorrect_noncarriers)/total_noncarriers < 0.3);
		if(report.num_correct_carriers + report.num_incorrect_noncarriers > 0.0){
			necessary_condition *= (((float) report.num_correct_carriers)/(report.num_correct_carriers + report.num_incorrect_noncarriers) > 0.7);
		}
		if(report.num_incorrect_carriers + report.num_correct_noncarriers > 0.0){
			necessary_condition *= (((float) report.num_correct_noncarriers)/(report.num_incorrect_carriers + report.num_correct_noncarriers) > 0.7);
		}
		if( necessary_condition && report.min > sum && report2.min > sum ){
			report.min         = sum;
			report.best_branch = node.label;
		}else{
			if( report.min > report2.min ){
				report.min         = report2.min;
				report.best_branch = report2.best_branch;
			}//else report is correct
		} 

		sum = report.num_correct_carriers + report.num_correct_noncarriers;

		necessary_condition = (((float) report.num_correct_carriers)/total_carriers < 0.3);
		necessary_condition *= (((float) report.num_correct_noncarriers)/total_noncarriers < 0.3);
		if(report.num_incorrect_carriers + report.num_correct_noncarriers > 0.0){
			necessary_condition *= (((float) report.num_incorrect_carriers)/(report.num_incorrect_carriers + report.num_correct_noncarriers) > 0.7);
		}
		if(report.num_correct_carriers + report.num_incorrect_noncarriers > 0.0){
			necessary_condition *= (((float) report.num_incorrect_noncarriers)/(report.num_correct_carriers + report.num_incorrect_noncarriers) > 0.7);
		}
		if( necessary_condition && report.flipped_min > sum && report2.flipped_min > sum ){
			report.flipped_min         = sum;
			report.best_flipped_branch = node.label;
		}else{
			if( report.flipped_min > report2.flipped_min ){
				report.flipped_min         = report2.flipped_min;
				report.best_flipped_branch = report2.best_flipped_branch;
			}//else report is correct
		}

	}else{

		if(sequences_carrying_mutations.member[node.label] == 1){
			report.num_correct_carriers      = 1;
			report.num_incorrect_carriers    = total_carriers - 1;
			report.num_correct_noncarriers   = total_noncarriers;
			report.num_incorrect_noncarriers = 0;

			if( report.num_incorrect_carriers/total_carriers < 0.3 ){
				report.min                     = report.num_incorrect_carriers;
				report.best_branch             = node.label;
			}else{
				report.min                     = std::numeric_limits<int>::max();
				report.best_branch             = -1;
			}
			if( report.num_correct_carriers/total_carriers < 0.3 && report.num_correct_noncarriers/total_noncarriers < 0.3 ){
				report.flipped_min             = report.num_correct_noncarriers + report.num_correct_carriers; 
				report.best_flipped_branch     = node.label;
			}else{
				report.flipped_min             = std::numeric_limits<int>::max();
				report.best_flipped_branch     = -1;
			} 
		}else{
			report.num_correct_carriers      = 0;
			report.num_incorrect_carriers    = total_carriers;
			report.num_correct_noncarriers   = total_noncarriers - 1;
			report.num_incorrect_noncarriers = 1;

			if( report.num_incorrect_carriers/total_carriers < 0.3 && report.num_incorrect_noncarriers/total_noncarriers < 0.3 ){
				report.min                     = report.num_incorrect_carriers + report.num_incorrect_noncarriers;
				report.best_branch             = node.label;
			}else{
				report.min                     = std::numeric_limits<int>::max();
				report.best_branch             = -1;
			}
			if( report.num_correct_noncarriers/total_noncarriers < 0.3 ){
				report.flipped_min             = report.num_correct_noncarriers;
				report.best_flipped_branch     = node.label;
			}else{
				report.flipped_min             = std::numeric_limits<int>::max();
				report.best_flipped_branch     = -1;
			} 
		}

	}

}

int 
MapMutation(Tree& tree, Leaves& sequences_carrying_mutations, Muts::iterator it_mut, bool use = false){

	(*it_mut).branch.clear();

	//if(sequences_carrying_mutations.num_leaves == 0 || sequences_carrying_mutations.num_leaves == N) return 1;
	if(sequences_carrying_mutations.num_leaves == 0) return 1;

	float min_value;
	int N = (tree.nodes.size() + 1.0)/2;
	float thr = (int) (0.03 * N);

	if(sequences_carrying_mutations.num_leaves == N){
		min_value = 0.0;
		(*it_mut).branch.resize(1);
		(*it_mut).flipped = false;
		(*it_mut).branch[0] = 2*N-2;
		tree.nodes[2*N-2].num_events += 1.0;
		return 1;
	}

	if(sequences_carrying_mutations.num_leaves == 0){
		min_value = 0.0;
		(*it_mut).branch.resize(0);
		(*it_mut).flipped = false;
		return 1;
	}

	PropagateStructGlobal report;
	PropagateMutation(*std::prev(tree.nodes.end(), 1), N, sequences_carrying_mutations, report);

	if(report.min == report.flipped_min && report.min <= thr){

		bool flag = true; //default flag
		if(flag){// not flipped
			min_value = report.min;
			(*it_mut).branch.resize(1);
			(*it_mut).branch[0] = report.best_branch; 
			(*it_mut).flipped = false;
			if(use) tree.nodes[report.best_branch].num_events += 1.0;
			assert((*it_mut).branch.size() == 1);
			return 1;
		}else{
			min_value = report.flipped_min;
			(*it_mut).branch.resize(1);
			(*it_mut).branch[0] = report.best_flipped_branch;
			(*it_mut).flipped = true;
			if(use) tree.nodes[report.best_flipped_branch].num_events += 1.0;
			assert((*it_mut).branch.size() == 1);
			return 2;
		}

	}else if( report.min <= report.flipped_min ){

		min_value = report.min;
		if( report.min <= thr ){
			(*it_mut).branch.resize(1);
			(*it_mut).branch[0] = report.best_branch; 
			(*it_mut).flipped = false;
			if(use) tree.nodes[report.best_branch].num_events += 1.0;
			assert((*it_mut).branch.size() == 1);
			return 1;
		}
		return 3;

	}else{

		min_value = report.flipped_min;
		if( report.flipped_min <= thr ){
			(*it_mut).branch.resize(1);
			(*it_mut).branch[0] = report.best_flipped_branch;
			(*it_mut).flipped = true;
			if(use) tree.nodes[report.best_flipped_branch].num_events += 1.0;
			assert((*it_mut).branch.size() == 1);
			return 2;
		}
		return 3;

	}

}


/////////////////////////////////

//convert to tree sequence (new set of nodes for each tree)
void
DumpAsTreeSequence(const std::string& filename_anc, const std::string& filename_mut, const std::string& filename_output){

	////////////////////////
	//read in anc file

	AncesTree anc;
	MarginalTree mtr, next_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	Muts::iterator it_mut, it_mut_current, it_mut_next, it_mut_tmp; //iterator for mut file
	float num_bases_tree_persists = 0.0, num_bases_next_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
	//The mut file is read once, file is closed after constructor is called.
	AncMutIterators ancmut(filename_anc, filename_mut);

	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	it_mut_current = it_mut;
	int N = (mtr.tree.nodes.size() + 1)/2.0, root = 2*N - 2, L = ancmut.NumSnps();
	int N_total = 2*N-1;
	Data data(N,L);
	std::vector<float> coordinates(2*data.N-1,0.0);

	Mutations mut;
	mut.Read(filename_mut);

	//........................................................................
	//Populate ts tables

	int ret;
	tsk_table_collection_t tables;
	ret = tsk_table_collection_init(&tables, 0);
	check_tsk_error(ret);

	tables.sequence_length = (*std::prev(ancmut.mut_end(),1)).pos + 1;
	for(int i = 0; i < N; i++){
		tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0, NULL, 0 , NULL, 0);
	}

	//population table

	//sites table
	char ancestral_allele[1];
	double pos, pos_begin, pos_end;
	std::vector<double> bps(L);
	int bps_index = 0;
	//tsk_site_table_add_row(&tables.sites, 1, ancestral_allele, sizeof(ancestral_allele), NULL, 0);
	for(; it_mut != ancmut.mut_end();){
		ancestral_allele[0] = (*it_mut).mutation_type[0];
		pos = (*it_mut).pos;
		int count = 0;

		it_mut_tmp = it_mut;
		while((*it_mut_tmp).pos == pos){
			it_mut_tmp++;
			count++;
			if(it_mut_tmp == ancmut.mut_end()) break;
		}
		assert(count > 0);

		if(count == 1){
			ret = tsk_site_table_add_row(&tables.sites, (*it_mut).pos, ancestral_allele, 1, NULL, 0);
			bps[bps_index] = (*it_mut).pos;
			bps_index++;
			it_mut++;
		}else{

			if(it_mut_tmp != ancmut.mut_end()){
				pos_end = ((*it_mut_tmp).pos + (*std::prev(it_mut_tmp)).pos)/2.0;
			}else{
				pos_end = (*std::prev(it_mut_tmp)).pos;
			}
			it_mut_tmp = it_mut;
			if(it_mut_tmp != it_mut_current){
				pos_begin = ((*it_mut_tmp).pos + (*std::prev(it_mut_tmp)).pos)/2.0;
			}else{
				pos_begin = pos;
			}
			int i = 0;
			while((*it_mut_tmp).pos == pos){
				ret = tsk_site_table_add_row(&tables.sites, ((i+1.0)/(count+1.0))*(pos_end - pos_begin) + pos_begin, ancestral_allele, 1, NULL, 0);
				bps[bps_index] = ((i+1.0)/(count+1.0))*(pos_end - pos_begin) + pos_begin;
				bps_index++;
				it_mut_tmp++;
				i++;
				if(it_mut_tmp == ancmut.mut_end()) break;
			}
			it_mut = it_mut_tmp;

		}

		//std::cerr << (*it_mut).pos << " " << count << " " << bps_index-1 << " " << bps[bps_index-1] << std::endl;
		check_tsk_error(ret);
	}
	assert(bps_index == L);

	if(ancmut.sample_ages.size() > 0){
		for(int i = 0; i < data.N; i++){
			ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, ancmut.sample_ages[i], TSK_NULL, i, NULL, 0);
			check_tsk_error(ret);
		}
	}else{
		for(int i = 0; i < data.N; i++){
			ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, i, NULL, 0);
			check_tsk_error(ret);
		}
	}

	///////////////////////////////////////////////////////////////////////////// 

	int snp, snp_end, tree_count = 0, node_rel, node, node_const, site_count = 0, count = 0;
	int node_count = data.N, edge_count = 0;

	std::vector<int> equivalent_branches_prev(N_total), equivalent_branches_next(N_total);
	float threshold_brancheq = 0.95;
	//float threshold_brancheq = 1.0;
	std::vector<std::vector<int>> potential_branches;
	//the number of leaves a branch needs to be equivalent
	potential_branches.resize(N);
	float threshold_inv = 1/(threshold_brancheq * threshold_brancheq);
	float N_float = N;
	for(int i = 1; i <= N; i++){
		potential_branches[i-1].push_back(i);
		//for branches with i+1 leaves, list the number of leaves a potential equivalent branch needs
		for(int j = i+1; j <= N; j++){
			if(threshold_inv >= j/(N_float-j) * ((N_float-i)/i) ){
				potential_branches[i-1].push_back(j);
				potential_branches[j-1].push_back(i);
			}
		}
	}

	num_bases_next_tree_persists = ancmut.NextTree(next_mtr, it_mut);
	it_mut_next = it_mut;

	////////////

	char derived_allele[1];
	char *branches;
	int prev_branch, next_branch;

	int metasize, size;
	char *meta;
	meta = (char *) malloc(1024);
	branches = (char *) malloc(1024);
	std::vector<int> SNPbegin(2*N-1,0.0),SNPend(2*N-1,0.0);
	while(num_bases_tree_persists >= 0.0){

		it_mut = it_mut_current;
		if(num_bases_next_tree_persists >= 0.0) anc.BranchAssociation(next_mtr.tree, mtr.tree, equivalent_branches_next, potential_branches, N, N_total, threshold_brancheq);

		mtr.tree.GetCoordinates(coordinates);
		for(int i = 0; i < mtr.tree.nodes.size()-1; i++){
			if(!(coordinates[(*mtr.tree.nodes[i].parent).label] - coordinates[i] > 0.0)){
				int parent = (*mtr.tree.nodes[i].parent).label, child = i;
				while(coordinates[parent] <= coordinates[child] + std::nextafter(coordinates[child], coordinates[child] + 1)){
					coordinates[parent] = coordinates[child] + std::nextafter(coordinates[child], coordinates[child] + 1);
					if(parent == root) break;
					child  = parent;
					parent = (*mtr.tree.nodes[parent].parent).label;
				}
			}
		}

		for(int i = 0; i < mtr.tree.nodes.size()-1; i++){  
			assert(coordinates[i] < coordinates[(*mtr.tree.nodes[i].parent).label]);
		}

		std::vector<int>::iterator it_snpbegin = SNPbegin.begin(), it_snpend = SNPend.begin();
		for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
			*it_snpbegin = ancmut.mut.info[(*it_node).SNP_begin].pos;
			if((*it_node).SNP_end < ancmut.mut.info.size()-1){
				*it_snpend   = ancmut.mut.info[(*it_node).SNP_end+1].pos;
			}else{
				*it_snpend   = ancmut.mut.info[(*it_node).SNP_end].pos;
			}
			it_snpbegin++;
			it_snpend++;
		}

		snp = mtr.pos;
		if(snp == 0){
			pos = 0;
		}else{
			pos = (bps[snp] + bps[snp-1])/2.0;
		}

		tree_count = (*it_mut).tree;
		node_const = tree_count * (data.N - 1);

		//Mutation table
		int l = snp;
		while((*it_mut).tree == tree_count){
			if((*it_mut).branch.size() == 1){
				node = *(*it_mut).branch.begin();
				if(node < N){
					derived_allele[0] = (*it_mut).mutation_type[2];
					ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, 0.5 * (coordinates[node] + coordinates[(*mtr.tree.nodes[node].parent).label]), derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}else{
					derived_allele[0] = (*it_mut).mutation_type[2];
					if(node < root){
						ret = tsk_mutation_table_add_row(&tables.mutations, l, node + node_const, TSK_NULL, 0.5*(coordinates[node] + coordinates[(*mtr.tree.nodes[node].parent).label]), derived_allele, 1, NULL, 0);
					}else{
						ret = tsk_mutation_table_add_row(&tables.mutations, l, node + node_const, TSK_NULL, coordinates[node], derived_allele, 1, NULL, 0);
					}
					check_tsk_error(ret);
				}
				site_count++;
			}

			l++;
			it_mut++; 
			if(l == L) break;
		}

		snp_end = l;
		if(snp_end < L){
			pos_end = (bps[snp_end-1] + bps[snp_end])/2.0;
		}else{
			pos_end = bps[L-1] + 1;
		}

		assert(pos != pos_end);
		assert(pos <= bps[snp]);
		assert(pos_end >= bps[snp]);

		//Node table
		int n = N;
		for(std::vector<float>::iterator it_coords = std::next(coordinates.begin(), data.N); it_coords != coordinates.end(); it_coords++){   

			if(count > 0){
				prev_branch = equivalent_branches_prev[n];
				if(prev_branch >= data.N) prev_branch += node_const - (data.N - 1);
			}else{
				prev_branch = -1;
			}
			if(num_bases_next_tree_persists >= 0.0){
				next_branch = equivalent_branches_next[n];
				if(next_branch >= data.N) next_branch += node_const + (data.N - 1);
			}else{
				next_branch = -1;
			}

			size = snprintf(NULL, 0,"%d",prev_branch) + snprintf(NULL, 0,"%d",next_branch) + 1;
			branches = (char *) realloc(branches, size);
			sprintf(branches, "%d %d", prev_branch, next_branch);

			ret = tsk_node_table_add_row(&tables.nodes, 0, *it_coords, TSK_NULL, TSK_NULL, branches, size);   
			check_tsk_error(ret);

			n++;
			node_count++;
		}

		std::vector<Node>::iterator it_node = std::next(mtr.tree.nodes.begin(), data.N);
		//Edge table
		for(it_node = mtr.tree.nodes.begin(); it_node != std::prev(mtr.tree.nodes.end(),1); it_node++){
			node = (*it_node).label;
			metasize = snprintf(NULL, 0,"%d",SNPbegin[node]) + snprintf(NULL, 0,"%d",SNPend[node]) + 1;
			meta = (char *) realloc(meta, metasize);
			sprintf(meta, "%d %d", SNPbegin[node], SNPend[node]);

			if(node >= data.N) node += node_const;

			if(0){
				if(count > 0){
					prev_branch = equivalent_branches_prev[(*it_node).label];
					if(prev_branch >= data.N) prev_branch += node_const - (data.N - 1);
				}else{
					prev_branch = -1;
				}
				if(num_bases_next_tree_persists >= 0.0){
					next_branch = equivalent_branches_next[(*it_node).label];
					if(next_branch >= data.N) next_branch += node_const + (data.N - 1);
				}else{
					next_branch = -1;
				}

				size = snprintf(NULL, 0,"%d",prev_branch) + snprintf(NULL, 0,"%d",next_branch) + 1;
				branches = (char *) realloc(branches, size);
				sprintf(branches, "%d %d", prev_branch, next_branch);
			}

			//ret = tsk_edge_table_add_row(&tables.edges, pos, pos_end, (*(*it_node).parent).label + node_const, node, branches, size);    
			ret = tsk_edge_table_add_row(&tables.edges, pos, pos_end, (*(*it_node).parent).label + node_const, node, meta, metasize);    
			check_tsk_error(ret);
			edge_count++;
		}

		//invert this vector
		std::fill(equivalent_branches_prev.begin(), equivalent_branches_prev.end(), -1);
		for(int i = 0; i < N_total; i++){
			if(equivalent_branches_next[i] != -1){
				equivalent_branches_prev[equivalent_branches_next[i]] = i;
			}
		}

		mtr = next_mtr;
		it_mut_current = it_mut_next;
		num_bases_tree_persists  = num_bases_next_tree_persists;
		if(num_bases_tree_persists >= 0){
			num_bases_next_tree_persists  = ancmut.NextTree(next_mtr, it_mut);
			it_mut_next = it_mut;
		}

		count++;

	} 
	//need to convert final tree


	ret = tsk_table_collection_sort(&tables, NULL, 0);
	check_tsk_error(ret);
	ret = tsk_table_collection_build_index(&tables, 0);
	check_tsk_error(ret);

	std::cerr << "Node count; edge count; tree count" << std::endl;
	std::cerr << node_count << " " << edge_count << " " << tree_count << std::endl;


	//////////////////////////

	// Write out the tree sequence
	ret = tsk_table_collection_dump(&tables, filename_output.c_str(), 0);        
	check_tsk_error(ret);
	tsk_table_collection_free(&tables); 

}

int
AddConstrainedNodeAgeLeastSquares(tsk_table_collection_t& tables, const int N, const double tolerance, const int iterations, const bool verbose){
	/* 
	 *  Use Dykstra's algorithm to solve the quadratic programming problem,
	 *
	 *    minimize \sum_i w_i (x_i - \hat{x}_i)^2
	 *    s.t. A x >= \epsilon
	 *
	 *  where `\hat{x}` are unconstrained node ages, `A` is a sparse matrix mapping node
	 *  ages onto branch lengths, `w` are weights, and `\epsilon` is the minimum allowed branch
	 *  length. In practice, we terminate if 0 < minimum_branch_length <= \epsilon.
	 *
	 *  NB: `tables` and `node_age` are modified _in place_.
	 */
	int max_iter = iterations; 
	double eps = 0.0;
	double tol = eps - tolerance; // stop if tol < min(branch_length) <= eps
	int ret = 1;

	assert (eps >= 0.0);
	assert (tol <= eps);

	std::vector<double> adj_len (tables.edges.num_rows, 0.0);

	tsk_size_t p, c;
	int iter = 0;
	double dual_gap, constraint;

	do {
		iter++;
		for (tsk_size_t e=0; e<tables.edges.num_rows; e++){
			p = tables.edges.parent[e];
			c = tables.edges.child[e];

			if (adj_len[e] > 0.0){ //constraint violated previously
				if (c >= N){
					tables.nodes.time[p] -= adj_len[e] * 0.5;
					tables.nodes.time[c] += adj_len[e] * 0.5;
				} else {
					tables.nodes.time[p] -= adj_len[e];
				}
			}

			adj_len[e] = tables.nodes.time[c] - tables.nodes.time[p] + eps;
			if (adj_len[e] > 0.0){ //constraint violated currently
				if (c >= N){ 
					tables.nodes.time[p] += adj_len[e] * 0.5;
					tables.nodes.time[c] -= adj_len[e] * 0.5;
				} else { 
					tables.nodes.time[p] += adj_len[e];
				}
			}
		}

		// residual for the dual problem
		dual_gap = std::numeric_limits<double>::infinity();
		for (tsk_size_t e=0; e<tables.edges.num_rows; e++){
			p = tables.edges.parent[e];
			c = tables.edges.child[e];
			constraint = tables.nodes.time[p] - tables.nodes.time[c];
			if (constraint < dual_gap){
				dual_gap = constraint;
			}
		}

		if (verbose && (iter % 10 == 0)) {
			std::cout << "\t[" << iter << "] min(branch length) = " << dual_gap << std::endl;
		}

		// success: min branch length is within convergence tol
		if (dual_gap > tol){
			if (verbose){
				std::cout << "Solution reached in " << iter << " iterations with " <<
					"min(branch length) = " << dual_gap << " > " << tol <<
					" ..." << std::endl;
			}
			ret = 0;
			break;
		}

		// failure: maximum iterations reached
		if (iter > max_iter){
			std::cerr << "Maximum number of iterations reached." << std::endl;
			ret = 1;
			break;
		}
	} while (true);

	return ret;
}

int
AddConstrainedNodeAgeBiased(tsk_table_collection_t& tables, const int N, const double tolerance){
	/* 
	 *  Constrain node age by moving parents backwards in time until they are
	 *  older than their children.
	 *
	 *  This relies on the tables being sorted, as the edges are ordered so that
	 *  edges with the same parent are together. `tsdate` uses the same strategy.
	 */

	double epsilon = tolerance;
	int ret = 0; 

	assert (epsilon >= 0.0);

	tsk_size_t parent = tables.edges.parent[0];
	tsk_size_t child = tables.edges.child[0];
	tsk_size_t last_parent = parent;
	double oldest_time = tables.nodes.time[child];
	for (tsk_size_t edge=0; edge<tables.edges.num_rows; edge++){
		parent = tables.edges.parent[edge];
		child = tables.edges.child[edge];
		if (parent == last_parent){
			if (oldest_time < tables.nodes.time[child]){
				oldest_time = tables.nodes.time[child];
			}
		} else {
			if (tables.nodes.time[last_parent] <= oldest_time){
				tables.nodes.time[last_parent] = oldest_time + epsilon;
			}
			oldest_time = tables.nodes.time[child];
		}
		last_parent = parent;
	}
	if (tables.nodes.time[last_parent] <= oldest_time){
		tables.nodes.time[last_parent] = oldest_time + epsilon;
	}

	return ret;
}

void
DumpAsCompressedTreeSequence(const std::string& filename_anc, const std::string& filename_mut, const std::string& filename_output, const double tolerance, const int iterations, const bool verbose){
	/* 
	 * Compress by combining equivalent branches across adjacent trees, and
	 * constraining node ages such that branch lengths are positive.
	 */

	MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	std::vector<Leaves> leaves, prev_leaves;
	Muts::iterator it_mut; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
	//The mut file is read once, file is closed after constructor is called.
	AncMutIterators ancmut(filename_anc, filename_mut);

	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	mtr.tree.FindAllLeaves(leaves); 
	int N = (mtr.tree.nodes.size() + 1)/2.0, root = 2*N - 2, L = ancmut.NumSnps(), T = ancmut.NumTrees();
	if (verbose){
		std::cout << "Read " << T << " marginal trees with " << N << " samples ..." << std::endl;
	}

	//........................................................................
	//Populate ts tables

	int ret;
	tsk_table_collection_t tables;
	ret = tsk_table_collection_init(&tables, 0);
	check_tsk_error(ret);

	// Add individual table
	tables.sequence_length = (*std::prev(ancmut.mut_end(),1)).pos + 1;
	for(int i = 0; i < N; i++){
		tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0, NULL, 0 , NULL, 0);
	}

	// Add population table
	//TODO

	// Add sites table
	char ancestral_allele[1];
	for(; it_mut != ancmut.mut_end(); it_mut++){
		ancestral_allele[0] = (*it_mut).mutation_type[0];
		ret = tsk_site_table_add_row(&tables.sites, (*it_mut).pos, ancestral_allele, 1, NULL, 0);
		check_tsk_error(ret);
	}

	//........................................................................
	//Iterate through ancmut

	it_mut = ancmut.mut_begin();
	std::vector<float> coordinates(2*N-1,0.0);
	int pos, snp, pos_end, snp_end, tree_count = 0, node, site_count = 0;

	// Numerator/denominator of average node age (weighted by tree span), stored in node table
	std::vector<double> node_span, node_age;
	node_span.reserve(N * T);
	node_age.reserve(N * T);
	double span, zero = 0.0, one = 1.0;

	// For each tree, keep a vector convert_nodes that maps marginal nodes to collapsed nodes
	int node_count = 0, edge_count = 0, root_count = 1;
	bool is_different = false;
	std::vector<int> convert_nodes(mtr.tree.nodes.size(), 0), convert_nodes_prev(mtr.tree.nodes.size(), 0);
	std::vector<int> update_backwards(2*N-1,0), update_forwards(2*N-1,0);

	std::vector<int>::iterator it_update_backwards = update_backwards.begin(), it_update_forwards = update_forwards.begin();
	std::vector<int>::iterator it_convert = convert_nodes.begin();
	for(; it_convert != std::next(convert_nodes.begin(),N); it_convert++){
		*it_convert = node_count;
		*it_update_backwards = node_count;
		*it_update_forwards  = node_count;

		if(ancmut.sample_ages.size() > 0){
			node_age.push_back((double)(ancmut.sample_ages[node_count]));
		} else {
			node_age.push_back(0.0);
		}
		node_span.push_back(1.0);
		ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, zero, TSK_NULL, TSK_NULL, (char*)(&zero), sizeof(double));   
		check_tsk_error(ret); 
		node_count++;
		it_update_forwards++;
		it_update_backwards++;
	}
	for(;it_convert != convert_nodes.end(); it_convert++){
		*it_convert = node_count;
		node_age.push_back(0.0);
		node_span.push_back(0.0);
		ret = tsk_node_table_add_row(&tables.nodes, 0, (double)(leaves[node_count].num_leaves-1), TSK_NULL, TSK_NULL, (char*)(&zero), sizeof(double));   
		check_tsk_error(ret); 
		node_count++;
	}

	double total_span = 0.0;
	char derived_allele[1];

	int metasize;
	char *meta;
	meta = (char *) malloc(1024);
	std::vector<float> prev_branch_persistence(2*N-1,0.0),branch_persistence(2*N-1,0.0);
	std::vector<int> SNPbegin(2*N-1,0.0),SNPend(2*N-1,0.0),prevSNPbegin(2*N-1,0.0),prevSNPend(2*N-1,0.0);
	while(num_bases_tree_persists >= 0.0){

		mtr.tree.GetCoordinates(coordinates);
		pos = (*it_mut).pos;
		if(mtr.pos == 0) pos = 0;

		prev_branch_persistence = branch_persistence;
		std::vector<float>::iterator it_branch_persistence = branch_persistence.begin();
		for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
			*it_branch_persistence = ancmut.mut.info[(*it_node).SNP_end].pos - ancmut.mut.info[(*it_node).SNP_begin].pos; //TODO:add mask file and add 0.5*dist to flanking SNPs
			it_branch_persistence++;
		}
		prevSNPbegin = SNPbegin;
		prevSNPend   = SNPend;
		std::vector<int>::iterator it_snpbegin = SNPbegin.begin(), it_snpend = SNPend.begin();
		for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
			*it_snpbegin = ancmut.mut.info[(*it_node).SNP_begin].pos;
			if((*it_node).SNP_end < ancmut.mut.info.size()-1){
				*it_snpend   = ancmut.mut.info[(*it_node).SNP_end+1].pos;
			}else{
				*it_snpend   = ancmut.mut.info[(*it_node).SNP_end].pos;
			}
			it_snpbegin++;
			it_snpend++;
		}


		for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node ++){
			(*it_node).SNP_begin = pos;
		}
		snp = mtr.pos;

		tree_count = (*it_mut).tree;

		if(tree_count > 0){
			//for each non-root node, check if its descendant set is identical to before
			//if not, update so that convert_nodes[i] = node_count
			std::fill(std::next(update_backwards.begin(),N), update_backwards.end(), 0);
			std::fill(std::next(update_forwards.begin(),N), update_forwards.end(), 0);
			std::fill(std::next(convert_nodes_prev.begin(),N), convert_nodes_prev.end(), 0);
			for(int n = 0; n < N; n++){
				mtr.tree.nodes[n].SNP_begin = prev_mtr.tree.nodes[n].SNP_begin;
			}
			//identify all nodes that are new (e.g. descendent set differs between prev and current tree)
			for(int n = N; n < 2*N-2; n++){
				if(1){
					is_different = true;
					if(leaves[n].num_leaves == prev_leaves[n].num_leaves){
						is_different = false;
						std::vector<int>::iterator it_leaves = leaves[n].member.begin(), it_prev_leaves = prev_leaves[n].member.begin();
						for(; it_leaves != leaves[n].member.end(); it_leaves++){
							if(*it_leaves != *it_prev_leaves){
								is_different = true;
								break;
							}
							it_prev_leaves++;
						}
					}
					if(is_different){ //still has a chance to be new
						if(1){ //for exact matching - should eventually be replaced by hash table
							for(int j = N; j < prev_leaves.size(); j++){
								if(leaves[n].num_leaves == prev_leaves[j].num_leaves){
									is_different = false;
									std::vector<int>::iterator it_prev_leaves = prev_leaves[j].member.begin();
									for(std::vector<int>::iterator it_leaves = leaves[n].member.begin(); it_leaves != leaves[n].member.end();){
										if(*it_leaves != *it_prev_leaves){
											is_different = true;
											break;
										}
										it_leaves++;
										it_prev_leaves++;
									}
									if(!is_different){ //found an identical node
										update_backwards[n]   = j; 
										update_forwards[j]    = n;
										convert_nodes_prev[n] = convert_nodes[j];
										mtr.tree.nodes[n].SNP_begin = prev_mtr.tree.nodes[j].SNP_begin;
										break;
									}
								}
							}
						}
					}else{
						update_backwards[n]   = n;
						update_forwards[n]    = n;
						convert_nodes_prev[n] = convert_nodes[n];
						mtr.tree.nodes[n].SNP_begin = prev_mtr.tree.nodes[n].SNP_begin;
					}
				}
			}

			// The root node will always have the same descendant set. However, if
			// the root is given the same ID across the entire ARG, it'll create a
			// strange constraint on branch lengths, because the TMRCA of all samples
			// will be forced to be constant across the sequence. To avoid this, give
			// a unique ID to the root whenever its children change.
			int root_left = (*mtr.tree.nodes[2*N-2].child_left).label,
					root_right = (*mtr.tree.nodes[2*N-2].child_right).label;
			is_different = update_backwards[root_left] == 0 || update_backwards[root_right] == 0;
			if (!is_different){
				update_backwards[2*N-2]   = 2*N-2;
				update_forwards[2*N-2]    = 2*N-2;
				convert_nodes_prev[2*N-2] = convert_nodes[2*N-2];
				mtr.tree.nodes[2*N-2].SNP_begin = prev_mtr.tree.nodes[2*N-2].SNP_begin;
			}
			root_count += int(is_different);

			// Update edge table
			for(int n = 0; n < 2*N-2; n++){
				int parent_prev = (*prev_mtr.tree.nodes[n].parent).label;
				int n_now       = update_forwards[n];
				int parent_now  = (*mtr.tree.nodes[n_now].parent).label;

				if(n < N){
					if( update_forwards[parent_prev] != parent_now ){
						//these edges don't exist anymore

						if(0){
							metasize = snprintf(NULL, 0,"%.2f",prev_branch_persistence[n]) + 1;
							meta = (char *) realloc(meta, metasize);
							sprintf(meta, "%.2f", prev_branch_persistence[n]);
						}
						metasize = snprintf(NULL, 0,"%d",prevSNPbegin[n]) + snprintf(NULL, 0,"%d",prevSNPend[n]) + 1;
						meta = (char *) realloc(meta, metasize);
						sprintf(meta, "%d %d", prevSNPbegin[n], prevSNPend[n]);

						if(0){
							if(prev_branch_persistence[n]+1000 < pos_end - prev_mtr.tree.nodes[n].SNP_begin){
								std::cerr << prev_branch_persistence[n] << " " << pos_end - prev_mtr.tree.nodes[n].SNP_begin << std::endl;
								std::cerr << "begin: " << SNPbegin[n] << " " << prev_mtr.tree.nodes[n].SNP_begin << std::endl;
								std::cerr << "end: " << SNPend[n] << " " << pos_end << std::endl;
							}
						}

						ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes[parent_prev], convert_nodes[n], meta, metasize);
						check_tsk_error(ret); 
						edge_count++;
						mtr.tree.nodes[n].SNP_begin = pos_end; 
					}
				}else if( n_now == 0 || update_forwards[parent_prev] != parent_now ){
					//these edges don't exist anymore

					if(0){
						metasize = snprintf(NULL, 0,"%.2f",prev_branch_persistence[n]) + 1;
						meta = (char *) realloc(meta, metasize);
						sprintf(meta, "%.2f", prev_branch_persistence[n]);
					}
					metasize = snprintf(NULL, 0,"%d",prevSNPbegin[n]) + snprintf(NULL, 0,"%d",prevSNPend[n]) + 1;
					meta = (char *) realloc(meta, metasize);
					sprintf(meta, "%d %d", prevSNPbegin[n], prevSNPend[n]);

					if(0){
						if(prev_branch_persistence[n]+1000 < pos_end - prev_mtr.tree.nodes[n].SNP_begin){
							std::cerr << prev_branch_persistence[n] << " " << pos_end - prev_mtr.tree.nodes[n].SNP_begin << std::endl;
							std::cerr << "begin: " << SNPbegin[n] << " " << prev_mtr.tree.nodes[n].SNP_begin << std::endl;
							std::cerr << "end: " << SNPend[n] << " " << pos_end << std::endl;
						}
					}

					ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes[parent_prev], convert_nodes[n], meta, metasize);
					if(n_now > 0) mtr.tree.nodes[n_now].SNP_begin = pos_end; 
					check_tsk_error(ret); 
					edge_count++; 
				}
			}

			// Update node table
			for(int n = N; n < 2*N-1; n++){
				if(update_backwards[n] == 0){
					convert_nodes[n] = node_count; //new node, so add to node table
					node_span.push_back(0.0);
					node_age.push_back(0.0);
					ret = tsk_node_table_add_row(&tables.nodes, 0, (double)(leaves[n].num_leaves-1), TSK_NULL, TSK_NULL, (char*)(&zero), sizeof(double));   
					mtr.tree.nodes[n].SNP_begin = pos; 
					node_count++;
				}else{
					convert_nodes[n] = convert_nodes_prev[n];
				}
			}          
		}

		// Update mutation table
		int l = snp;
		while((*it_mut).tree == tree_count){
			if((*it_mut).branch.size() == 1){
				node = *(*it_mut).branch.begin();
				if(node < N){
					derived_allele[0] = (*it_mut).mutation_type[2];
					ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, TSK_UNKNOWN_TIME, derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}else{
					derived_allele[0] = (*it_mut).mutation_type[2];
					ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes[node], TSK_NULL, TSK_UNKNOWN_TIME, derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}
				site_count++;
			}

			l++;
			it_mut++; 
			if(l == L) break;
		}
		snp_end = l;
		if(snp_end < L){
			pos_end = (*it_mut).pos;
		} else {
			pos_end = (*std::prev(ancmut.mut_end(),1)).pos + 1;
		}

		// Update average node age and span 
		// (the latter stored in node metadata; it is a char array so cast to double for arithmetic)
		span = pos_end - pos;
		for(int n = N; n < 2*N-1; n++){
			node_age[convert_nodes[n]] += span * (double)(coordinates[n]); //TODO overflow prone
			node_span[convert_nodes[n]] += span;
		}
		total_span += span;

		// Progress
		int chunk_size = int(T * 0.05);
		int percent_complete = int(100 * double(tree_count) / double(T));
		if (verbose && tree_count % chunk_size == 0)
		{
			std::cout << "\tProcessed " << node_count << "/" << edge_count << " nodes/edges from " << 
				tree_count << " trees (" << percent_complete << "%)" << std::endl;
		}

		// Load next tree
		prev_mtr                = mtr;
		prev_leaves             = leaves;
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		mtr.tree.FindAllLeaves(leaves);
	} 

	// For last tree, dump all edges
	tree_count++;
	for(int n = 0; n < 2*N-2; n++){        
		int parent_prev = (*prev_mtr.tree.nodes[n].parent).label;

		if(0){
			metasize = snprintf(NULL, 0,"%.2f",branch_persistence[n]) + 1;
			meta = (char *) realloc(meta, metasize);
			sprintf(meta, "%.2f", branch_persistence[n]);
		}
		metasize = snprintf(NULL, 0,"%d",SNPbegin[n]) + snprintf(NULL, 0,"%d",prevSNPend[n]) + 1;
		meta = (char *) realloc(meta, metasize);
		sprintf(meta, "%d %d", SNPbegin[n], SNPend[n]);

		if(0){
			if(branch_persistence[n]+1000 < pos_end - prev_mtr.tree.nodes[n].SNP_begin){
				std::cerr << branch_persistence[n] << " " << pos_end - prev_mtr.tree.nodes[n].SNP_begin << std::endl;
			}
		}

		ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes[parent_prev], convert_nodes[n], meta, metasize);
		check_tsk_error(ret); 
		edge_count++;
	}

	if (verbose){
		std::cout << "\tNodes: " << node_count << "\n" <<
			"\tEdges: " << edge_count << "\n" <<
			"\tTrees: " << tree_count << "\n" <<
			"\tRoots: " << root_count << std::endl;
	}

	// Copy node ages into metadata
	assert (node_age.size() == tables.nodes.num_rows);
	assert (node_span.size() == tables.nodes.num_rows);
	if (verbose){
		std::cout << "Storing unconstrained node age as double precision in node metadata ..." << std::endl;
	}
	for (tsk_size_t i=0; i<tables.nodes.num_rows; i++){
		node_age[i] /= node_span[i];
	}
	std::memcpy(tables.nodes.metadata, node_age.data(), sizeof(double)*tables.nodes.num_rows);
	node_age.clear();
	node_age.shrink_to_fit();
	node_span.clear();
	node_span.shrink_to_fit();

	// Sort table (needed for edge order)
	if (verbose){ 
		std::cout << "Sorting tree sequence table collection ..." << std::endl;
	}
	ret = tsk_table_collection_sort(&tables, NULL, 0);
	check_tsk_error(ret);
	ret = tsk_table_collection_build_index(&tables, 0);
	check_tsk_error(ret);

	// Find closest node ages (in least-squares sense) that result in
	// positive branch lengths.
	if (verbose){ 
		std::cout << "Constraining node ages ..." << std::endl;
	}
	std::memcpy(tables.nodes.time, tables.nodes.metadata, sizeof(double)*tables.nodes.num_rows);
	if (iterations > 0){
		ret = AddConstrainedNodeAgeLeastSquares(tables, N, tolerance, iterations, verbose);
	}
	ret = AddConstrainedNodeAgeBiased(tables, N, tolerance); // force constraint by pushing node ages backwards in time

	// Sort table
	if (verbose){ 
		std::cout << "Resorting tree sequence table collection ..." << std::endl;
	}
	ret = tsk_table_collection_sort(&tables, NULL, 0);
	check_tsk_error(ret);
	ret = tsk_table_collection_build_index(&tables, 0);
	check_tsk_error(ret);

	//////////////////////////

	// Write out the tree sequence
	if (verbose){ 
		std::cout << "Writing tree sequence ..." << std::endl;
	}
	ret = tsk_table_collection_dump(&tables, filename_output.c_str(), 0);        
	check_tsk_error(ret);
	tsk_table_collection_free(&tables); 

	if (verbose){ 
		std::cout << "Finished." << std::endl;
	}
}

void
FindIdenticalNodes(Tree& prev_tr, Tree& tr, std::vector<Leaves>& leaves, std::vector<Leaves>& prev_leaves, std::vector<int>& update_forwards, std::vector<int>& update_backwards, std::vector<int>& prev_rewire, std::vector<int>& rewire, std::vector<int>& convert_nodes, std::vector<int>& convert_nodes_prev, int N){

	//nodes:
	//for each node, check if its descendant set is identical to before
	//if no, update convert_nodes[i] = node_count;
	std::fill(std::next(update_backwards.begin(),N), update_backwards.end(), 0);
	std::fill(std::next(update_forwards.begin(),N), update_forwards.end(), 0);

	for(int n = 0; n < N; n++){
		tr.nodes[n].SNP_begin = prev_tr.nodes[n].SNP_begin;
	}

	//identify all nodes that are new in mtr compared to prev_mtr
	//TODO: use hash table
	bool is_different = true;
	for(int n = N; n < 2*N-1; n++){
		is_different = true;
		if(leaves[n].num_leaves == prev_leaves[n].num_leaves){
			is_different = false;
			std::vector<int>::iterator it_leaves = leaves[n].member.begin(), it_prev_leaves = prev_leaves[n].member.begin();
			for(; it_leaves != leaves[n].member.end(); it_leaves++){
				if(*it_leaves != *it_prev_leaves){
					is_different = true;
					break;
				}
				it_prev_leaves++;
			}
		}
		if(is_different){ //still has a chance to be new

			if(1){ //for exact matching - should eventually be replaced by hash table
				for(int j = N; j < prev_leaves.size(); j++){
					if(leaves[n].num_leaves == prev_leaves[j].num_leaves){

						is_different = false;
						std::vector<int>::iterator it_prev_leaves = prev_leaves[j].member.begin();
						for(std::vector<int>::iterator it_leaves = leaves[n].member.begin(); it_leaves != leaves[n].member.end();){
							if(*it_leaves != *it_prev_leaves){
								is_different = true;
								break;
							}
							it_leaves++;
							it_prev_leaves++;
						}

						if(convert_nodes_prev[j] != 0) assert(prev_rewire[j] == j);
						if(!is_different){ //found an identical node
							update_backwards[n]   = j; 
							update_forwards[j]    = n;
							tr.nodes[n].SNP_begin = prev_tr.nodes[j].SNP_begin;
							if(prev_rewire[j] == j){
								assert(convert_nodes_prev[j] != 0);
								convert_nodes[n] = convert_nodes_prev[j];
								rewire[n] = n;
							}else{
								//new node
								prev_rewire[j]        = j;
								rewire[n]             = n;
								assert(convert_nodes_prev[j] == 0);
							}
							break;
						} 

					}
				}
			}

		}else{

			update_backwards[n]   = n;
			update_forwards[n]    = n;
			tr.nodes[n].SNP_begin = prev_tr.nodes[n].SNP_begin;
			if(prev_rewire[n] == n){
				assert(convert_nodes_prev[n] != 0);
				convert_nodes[n] = convert_nodes_prev[n];
				rewire[n] = n;
			}else{
				//new node
				prev_rewire[n]        = n;
				rewire[n]             = n;
				assert(convert_nodes_prev[n] == 0);
			}

		}
	}


}

//removes branches with no mutation mapped to it (TODO: this is not yet optimal in terms of compression)
void
DumpAsTreeSequenceWithPolytomies(const std::string& filename_anc, const std::string& filename_mut, const std::string& filename_output){

	MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	std::vector<Leaves> leaves, prev_leaves;
	Muts::iterator it_mut, it_mut_prev_begin, it_mut_begin; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
	//The mut file is read once, file is closed after constructor is called.
	AncMutIterators ancmut(filename_anc, filename_mut);

	num_bases_tree_persists = ancmut.NextTree(prev_mtr, it_mut_prev_begin);
	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut_begin);
	prev_mtr.tree.FindAllLeaves(prev_leaves);
	mtr.tree.FindAllLeaves(leaves);
	int N = (mtr.tree.nodes.size() + 1)/2.0, root = 2*N - 2, L = ancmut.NumSnps();

	//........................................................................
	//Populate ts tables

	int ret;
	tsk_table_collection_t tables;
	ret = tsk_table_collection_init(&tables, 0);
	check_tsk_error(ret);

	tables.sequence_length = (*std::prev(ancmut.mut_end(),1)).pos + 1;
	for(int i = 0; i < N; i++){
		tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0, NULL, 0 , NULL, 0);
	}

	//population table

	//sites table
	char ancestral_allele[1];
	//tsk_site_table_add_row(&tables.sites, 1, ancestral_allele, sizeof(ancestral_allele), NULL, 0);
	for(it_mut = ancmut.mut_begin(); it_mut != ancmut.mut_end(); it_mut++){
		ancestral_allele[0] = (*it_mut).mutation_type[0];
		ret = tsk_site_table_add_row(&tables.sites, (*it_mut).pos, ancestral_allele, 1, NULL, 0);
		check_tsk_error(ret);
	}

	/////////////////////////////////////////////

	std::vector<int> update_backwards(2*N-1,0), update_forwards(2*N-1,0);
	std::vector<int> rewire(2*N-1,0), prev_rewire(2*N-1,0);
	//for each tree, keep a vector convert_nodes that maps nodes to ts nodes
	std::vector<int> convert_nodes(2*N-1, 0), convert_nodes_prev(2*N-1, 0);

	int node_count = N, edge_count = 0;
	prev_rewire[root]        = root;
	convert_nodes_prev[root] = node_count;
	rewire[root]             = root;
	convert_nodes[root]      = node_count; 
	node_count++; 
	for(std::vector<Node>::iterator it_node = prev_mtr.tree.nodes.begin(); it_node != prev_mtr.tree.nodes.end(); it_node++){
		(*it_node).SNP_begin = 0;
	}
	for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
		(*it_node).SNP_begin = (*it_mut_begin).pos;
	}

	for(int i = 0; i < N; i++){
		prev_rewire[i]        = i;
		rewire[i]             = i;
		convert_nodes_prev[i] = i;
		convert_nodes[i]      = i;
		update_forwards[i]    = i;
		update_backwards[i]   = i;
		ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, prev_leaves[i].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
		check_tsk_error(ret);
	}

	ret = tsk_node_table_add_row(&tables.nodes, 0, prev_leaves[root].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
	check_tsk_error(ret);
	for(int i = N; i < 2*N-2; i++){
		if(prev_mtr.tree.nodes[i].num_events >= 1.0){
			prev_rewire[i]        = i;
			convert_nodes_prev[i] = node_count;
			ret = tsk_node_table_add_row(&tables.nodes, 0, prev_leaves[i].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);   
			check_tsk_error(ret);
			node_count++;
		}
		if(mtr.tree.nodes[i].num_events >= 1.0){
			rewire[i] = i;
			//don't know if this is a new node yet
		}
	}

	//find identical nodes
	FindIdenticalNodes(prev_mtr.tree, mtr.tree, leaves, prev_leaves, update_forwards, update_backwards, prev_rewire, rewire, convert_nodes, convert_nodes_prev, N);

	for(int i = N; i < 2*N-2; i++){
		if(prev_mtr.tree.nodes[i].num_events < 1.0){
			int k = i;
			while(prev_rewire[k] != k){
				k         = (*prev_mtr.tree.nodes[k].parent).label; 
				prev_rewire[i] = k;
				if(k == root){
					prev_rewire[i] = root;
					break;
				}
			}
		}else{
			assert(prev_rewire[i] == i);
		}
	}

	//add new nodes
	for(int n = N; n < 2*N-2; n++){
		if(rewire[n] == n && update_backwards[n] == 0){
			//new node
			convert_nodes[n] = node_count;
			ret = tsk_node_table_add_row(&tables.nodes, 0, leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);   
			check_tsk_error(ret);
			mtr.tree.nodes[n].SNP_begin = (*it_mut_begin).pos;
			node_count++;
		}
		if(prev_rewire[n] == n && convert_nodes_prev[n] == 0){
			//new node
			convert_nodes_prev[n] = node_count;
			prev_mtr.tree.nodes[n].SNP_begin = 0;
			assert(update_forwards[n] != 0);
			if(update_forwards[n] != 0){
				assert(convert_nodes[update_forwards[n]] == 0);
				convert_nodes[update_forwards[n]] = node_count;
				mtr.tree.nodes[update_forwards[n]].SNP_begin = 0;
			}
			ret = tsk_node_table_add_row(&tables.nodes, 0, prev_leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
			check_tsk_error(ret);
			node_count++;
		}
		if(prev_rewire[n] == n){
			assert(convert_nodes_prev[n] > 0);
		}
		if(convert_nodes_prev[n] != 0){
			assert(prev_rewire[n] == n);
		}
	}

	//........................................................................
	//Iterate through ancmut
	//
	it_mut = ancmut.mut_begin();
	int pos, snp, pos_end, snp_end, tree_count = 0, node, site_count = 0;
	bool is_different = false;
	char derived_allele[1];
	while(num_bases_tree_persists >= 0.0){

		//mtr.tree.GetCoordinates(coordinates);
		pos = (*it_mut).pos;
		if(prev_mtr.pos == 0) pos = 0;
		snp = prev_mtr.pos;
		tree_count = (*it_mut).tree;

		//Mutation table
		int l = snp;
		while((*it_mut).tree == tree_count){
			if((*it_mut).branch.size() == 1 && (*it_mut).flipped == 0){
				node = *(*it_mut).branch.begin();
				if(node < N){
					derived_allele[0] = (*it_mut).mutation_type[2];
					ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, 0, derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}else{
					derived_allele[0] = (*it_mut).mutation_type[2];
					assert(prev_rewire[node] == node);
					assert(convert_nodes_prev[node] != 0);
					ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes_prev[node], TSK_NULL, 0, derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}
				site_count++;
			}

			l++;
			it_mut++; 
			if(l == L) break;
		}
		snp_end = l;
		if(snp_end < L){
			pos_end = (*it_mut).pos;
		}else{
			pos_end = (*std::prev(ancmut.mut_end(),1)).pos + 1;
		}

		//write edges

		//node is new if update[n] == 0 and rewire[n] == n
		//if rewire[n] != n, this node is removed
		//if update[n] > 0 and rewire[n] == n, then it is identical to a node in the prev tree
		for(int n = 0; n < 2*N-2; n++){
			int parent_prev = prev_rewire[(*prev_mtr.tree.nodes[n].parent).label];
			int n_now       = update_forwards[n];
			int parent_now  = rewire[(*mtr.tree.nodes[n_now].parent).label]; 
			//rewire might still change, however if it was identical to parent_prev, this would be already reflected here

			if(n < N){
				if( update_forwards[parent_prev] == 0 || update_forwards[parent_prev] != parent_now ){
					//these edges don't exist anymore
					if(n > 0) assert(convert_nodes_prev[n] != 0);
					assert(convert_nodes_prev[parent_prev] != 0);
					ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n], NULL, 0);
					check_tsk_error(ret); 
					edge_count++;
					//prev_mtr.tree.nodes[n].SNP_begin = pos_end;
					mtr.tree.nodes[n].SNP_begin = pos_end; 
				}
			}else if( prev_rewire[n] == n && (n_now == 0 || update_forwards[parent_prev] == 0 || parent_now == 0 || update_forwards[parent_prev] != parent_now) ){
				//these edges don't exist anymore
				assert(convert_nodes_prev[n] != 0);
				assert(convert_nodes_prev[parent_prev] != 0);
				ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n], NULL, 0);
				check_tsk_error(ret);
				if(n_now > 0) mtr.tree.nodes[n_now].SNP_begin = pos_end; 
				edge_count++; 
			}
		}

		//load next tree
		prev_mtr                = mtr;
		prev_leaves             = leaves;
		prev_rewire             = rewire;
		convert_nodes_prev      = convert_nodes;
		it_mut_prev_begin       = it_mut_begin;
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

		//skip last part if I have just read the last tree
		if(num_bases_tree_persists >= 0){    

			it_mut_begin            = it_mut;
			it_mut                  = it_mut_prev_begin;
			std::fill(std::next(rewire.begin(),N), std::prev(rewire.end(),1), 0);
			std::fill(std::next(convert_nodes.begin(),N), std::prev(convert_nodes.end(),1), 0);
			mtr.tree.FindAllLeaves(leaves);
			for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node ++){
				(*it_node).SNP_begin = (*it_mut_begin).pos;
			}
			for(int i = 0; i < 2*N-2; i++){
				if(mtr.tree.nodes[i].num_events >= 1.0){
					rewire[i] = i;
				}
			}

			//identify equivalent branches,
			//populate update_forwards and update_backwards
			//set prev_rewire[i] = i and rewire[k] = k if these nodes should exist
			//transfer relevant nodes from convert_nodes_prev to convert_nodes 
			FindIdenticalNodes(prev_mtr.tree, mtr.tree, leaves, prev_leaves, update_forwards, update_backwards, prev_rewire, rewire, convert_nodes, convert_nodes_prev, N);
			for(int i = N; i < 2*N-2; i++){
				if(prev_mtr.tree.nodes[i].num_events < 1.0){
					int k = i;
					while(prev_rewire[k] != k){
						k         = (*prev_mtr.tree.nodes[k].parent).label; 
						prev_rewire[i] = k;
						if(k == root){
							prev_rewire[i] = root;
							break;
						}
					}
				}else{
					assert(prev_rewire[i] == i);
				}
			}

			//write new nodes
			for(int n = N; n < 2*N-2; n++){
				if(update_backwards[n] == 0 && rewire[n] == n){
					convert_nodes[n] = node_count;
					//new node, so add to node table
					//need to think how I can replace num_leaves by actual age
					ret = tsk_node_table_add_row(&tables.nodes, 0, leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
					check_tsk_error(ret); 
					mtr.tree.nodes[n].SNP_begin = (*it_mut_begin).pos; 
					node_count++;
				}
				if(prev_rewire[n] == n && convert_nodes_prev[n] == 0){
					convert_nodes_prev[n] = node_count;
					//if(update_forwards[n] == 0) std::cerr << n << std::endl;
					prev_mtr.tree.nodes[n].SNP_begin = (*it_mut_prev_begin).pos;
					assert(update_forwards[n] != 0);
					if(update_forwards[n] != 0){
						convert_nodes[update_forwards[n]] = node_count;
						mtr.tree.nodes[update_forwards[n]].SNP_begin = prev_mtr.tree.nodes[n].SNP_begin;
					}
					ret = tsk_node_table_add_row(&tables.nodes, 0, prev_leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0); 
					check_tsk_error(ret); 
					node_count++;
				}
				if(prev_rewire[n] == n){
					assert(convert_nodes_prev[n] > 0);
				}
				if(convert_nodes_prev[n] != 0){
					assert(prev_rewire[n] == n);
				}
			}      

		}   

	} 


	for(int i = N; i < 2*N-2; i++){
		if(prev_mtr.tree.nodes[i].num_events < 1.0){
			int k = i;
			while(prev_rewire[k] != k){
				k         = (*prev_mtr.tree.nodes[k].parent).label; 
				prev_rewire[i] = k;
				if(k == root){
					prev_rewire[i] = root;
					break;
				}
			}
		}else{
			assert(prev_rewire[i] == i);
		}
	}

	//last tree
	pos = (*it_mut_prev_begin).pos;
	snp = prev_mtr.pos;
	tree_count = (*it_mut).tree;
	//Mutation table
	int l = snp;
	while((*it_mut).tree == tree_count){
		if((*it_mut).branch.size() == 1 && (*it_mut).flipped == 0){
			node = *(*it_mut).branch.begin();
			if(node < N){
				derived_allele[0] = (*it_mut).mutation_type[2];
				ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, 0, derived_allele, 1, NULL, 0);
				check_tsk_error(ret);
			}else{
				derived_allele[0] = (*it_mut).mutation_type[2];
				assert(prev_rewire[node] == node);
				assert(convert_nodes_prev[node] != 0);
				ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes_prev[node], TSK_NULL, 0, derived_allele, 1, NULL, 0);
				check_tsk_error(ret);
			}
			site_count++;
		}

		l++;
		it_mut++; 
		if(l == L) break;
	}
	snp_end = l;
	if(snp_end < L){
		pos_end = (*it_mut).pos;
	}else{
		pos_end = (*std::prev(ancmut.mut_end(),1)).pos + 1;
	}

	//for last tree need to dump all edges
	for(int n = 0; n < 2*N-2; n++){
		//Edge table
		//these edges don't exist anymore
		if(rewire[n] == n){

			if(n > 0) assert(convert_nodes_prev[n] != 0);
			int parent_prev = prev_rewire[(*prev_mtr.tree.nodes[n].parent).label];
			assert(convert_nodes_prev[parent_prev] != 0);
			ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n], NULL, 0);   
			check_tsk_error(ret);
			edge_count++;
		}
	}

	std::cerr << "Node count; edge count; tree count" << std::endl;
	std::cerr << node_count << " " << edge_count << " " << tree_count << std::endl;

	tsk_table_collection_sort(&tables, NULL, 0);
	check_tsk_error(ret);
	tsk_table_collection_build_index(&tables, 0);
	check_tsk_error(ret);
	//////////////////////////

	// Write out the tree sequence
	ret = tsk_table_collection_dump(&tables, filename_output.c_str(), 0);        
	check_tsk_error(ret);
	tsk_table_collection_free(&tables); 

}

//removes branches with no mutation mapped to it, where mutations are remapped so that data can be recovered exactly (TODO: this is not yet optimal in terms of compression)
void
DumpAsTreeSequenceWithPolytomies(const std::string& filename_anc, const std::string& filename_mut, const std::string& filename_haps, const std::string& filename_sample, const std::string& filename_output){

	haps m_hap(filename_haps.c_str(), filename_sample.c_str());
	Data data(m_hap.GetN(), m_hap.GetL());

	Leaves sequences_carrying_mutation;
	sequences_carrying_mutation.member.resize(data.N);
	std::vector<char> sequence(data.N);
	int bp;


	MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	std::vector<Leaves> leaves, prev_leaves;
	Muts::iterator it_mut, it_mut_prev_begin, it_mut_begin; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
	//The mut file is read once, file is closed after constructor is called.
	AncMutIterators ancmut(filename_anc, filename_mut);

	num_bases_tree_persists = ancmut.NextTree(prev_mtr, it_mut_prev_begin);
	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut_begin);
	prev_mtr.tree.FindAllLeaves(prev_leaves);
	mtr.tree.FindAllLeaves(leaves);
	int N = (mtr.tree.nodes.size() + 1)/2.0, root = 2*N - 2, L = ancmut.NumSnps();

	//........................................................................
	//Populate ts tables

	int ret;
	tsk_table_collection_t tables;
	ret = tsk_table_collection_init(&tables, 0);
	check_tsk_error(ret);

	tables.sequence_length = (*std::prev(ancmut.mut_end(),1)).pos + 1;
	for(int i = 0; i < N; i++){
		tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0, NULL, 0 , NULL, 0);
	}

	//population table

	//sites table
	char ancestral_allele[1];
	//tsk_site_table_add_row(&tables.sites, 1, ancestral_allele, sizeof(ancestral_allele), NULL, 0);
	for(it_mut = ancmut.mut_begin(); it_mut != ancmut.mut_end(); it_mut++){
		ancestral_allele[0] = (*it_mut).mutation_type[0];
		ret = tsk_site_table_add_row(&tables.sites, (*it_mut).pos, ancestral_allele, 1, NULL, 0);
		check_tsk_error(ret);
	}


	///////////////////////////////////
	//remap
	for(std::vector<Node>::iterator it_node = prev_mtr.tree.nodes.begin(); it_node != prev_mtr.tree.nodes.end(); it_node++){
		(*it_node).num_events = 0.0;
	}
	//remap mutations
	it_mut = ancmut.mut_begin();
	while((*it_mut).tree == ancmut.get_treecount()-1){
		m_hap.ReadSNP(sequence, bp);
		if(bp != (*it_mut).pos){
			std::cerr << "Error: haps file and anc/mut files don't contain same set of SNPs." << std::endl;
			exit(1);
		}

		sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
		for(int i = 0; i < data.N; i++){
			if(sequence[i] == '1'){
				sequences_carrying_mutation.member[i] = 1;
				sequences_carrying_mutation.num_leaves++;
			}else{
				sequences_carrying_mutation.member[i] = 0;
			}
		}

		if(sequences_carrying_mutation.num_leaves > 0 && sequences_carrying_mutation.num_leaves < N){
			int num_b = MapMutationExact(prev_mtr.tree, sequences_carrying_mutation, it_mut);
		}
		it_mut++;
	}
	while((*it_mut).tree == ancmut.get_treecount()){
		m_hap.ReadSNP(sequence, bp);
		if(bp != (*it_mut).pos){
			std::cerr << "Error: haps file and anc/mut files don't contain same set of SNPs." << std::endl;
			exit(1);
		}

		sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
		for(int i = 0; i < data.N; i++){
			if(sequence[i] == '1'){
				sequences_carrying_mutation.member[i] = 1;
				sequences_carrying_mutation.num_leaves++;
			}else{
				sequences_carrying_mutation.member[i] = 0;
			}
		}

		if(sequences_carrying_mutation.num_leaves > 0 && sequences_carrying_mutation.num_leaves < N){
			MapMutationExact(mtr.tree, sequences_carrying_mutation, it_mut);
		}
		it_mut++;
	}
	it_mut = ancmut.mut_begin();

	/////////////////////////////////////////////

	std::vector<int> update_backwards(2*N-1,0), update_forwards(2*N-1,0);
	std::vector<int> rewire(2*N-1,0), prev_rewire(2*N-1,0);
	//for each tree, keep a vector convert_nodes that maps nodes to ts nodes
	std::vector<int> convert_nodes(2*N-1, 0), convert_nodes_prev(2*N-1, 0);

	int node_count = N, edge_count = 0;
	prev_rewire[root] = root;
	convert_nodes_prev[root] = node_count;
	rewire[root] = root;
	convert_nodes[root] = node_count; 
	node_count++; 
	for(std::vector<Node>::iterator it_node = prev_mtr.tree.nodes.begin(); it_node != prev_mtr.tree.nodes.end(); it_node++){
		(*it_node).SNP_begin = 0;
	}
	for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
		(*it_node).SNP_begin = (*it_mut_begin).pos;
	}

	for(int i = 0; i < N; i++){
		prev_rewire[i]        = i;
		rewire[i]             = i;
		convert_nodes_prev[i] = i;
		convert_nodes[i]      = i;
		update_forwards[i]    = i;
		update_backwards[i]   = i;
		ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, prev_leaves[i].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
		check_tsk_error(ret);
	}

	ret = tsk_node_table_add_row(&tables.nodes, 0, prev_leaves[root].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
	check_tsk_error(ret);
	for(int i = N; i < 2*N-2; i++){
		if(prev_mtr.tree.nodes[i].num_events >= 1.0){
			prev_rewire[i] = i;
			convert_nodes_prev[i] = node_count;
			ret = tsk_node_table_add_row(&tables.nodes, 0, prev_leaves[i].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);   
			check_tsk_error(ret);
			node_count++;
		}
		if(mtr.tree.nodes[i].num_events >= 1.0){
			rewire[i] = i;
			//don't know if this is a new node yet
		}
	}

	//find identical nodes
	FindIdenticalNodes(prev_mtr.tree, mtr.tree, leaves, prev_leaves, update_forwards, update_backwards, prev_rewire, rewire, convert_nodes, convert_nodes_prev, N);

	for(int i = N; i < 2*N-2; i++){
		if(prev_mtr.tree.nodes[i].num_events < 1.0){
			int k = i;
			while(prev_rewire[k] != k){
				k         = (*prev_mtr.tree.nodes[k].parent).label; 
				prev_rewire[i] = k;
				if(k == root){
					prev_rewire[i] = root;
					break;
				}
			}
		}else{
			assert(prev_rewire[i] == i);
		}
	}

	//add new nodes
	for(int n = N; n < 2*N-2; n++){
		if(rewire[n] == n && update_backwards[n] == 0){
			//new node
			convert_nodes[n] = node_count;
			ret = tsk_node_table_add_row(&tables.nodes, 0, leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);   
			check_tsk_error(ret);
			mtr.tree.nodes[n].SNP_begin = (*it_mut_begin).pos;
			node_count++;
		}
		if(prev_rewire[n] == n && convert_nodes_prev[n] == 0){
			convert_nodes_prev[n] = node_count;
			prev_mtr.tree.nodes[n].SNP_begin = 0;
			assert(update_forwards[n] != 0);
			if(update_forwards[n] != 0){
				convert_nodes[update_forwards[n]] = node_count;
				mtr.tree.nodes[update_forwards[n]].SNP_begin = prev_mtr.tree.nodes[n].SNP_begin;
			}
			ret = tsk_node_table_add_row(&tables.nodes, 0, prev_leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
			check_tsk_error(ret);
			node_count++;
		}
		if(prev_rewire[n] == n){
			assert(convert_nodes_prev[n] > 0);
		}
		if(convert_nodes_prev[n] != 0){
			assert(prev_rewire[n] == n);
		}
	}

	//........................................................................
	//Iterate through ancmut
	//
	it_mut = ancmut.mut_begin();
	int pos, snp, pos_end, snp_end, tree_count = 0, node, site_count = 0;
	bool is_different = false;
	char derived_allele[1];
	while(num_bases_tree_persists >= 0.0){

		//mtr.tree.GetCoordinates(coordinates);
		pos = (*it_mut).pos;
		if(prev_mtr.pos == 0) pos = 0;
		snp = prev_mtr.pos;
		tree_count = (*it_mut).tree;

		//Mutation table
		int l = snp;
		while((*it_mut).tree == tree_count){
			for(std::vector<int>::iterator it_branch = (*it_mut).branch.begin(); it_branch != (*it_mut).branch.end(); it_branch++){ 
				node = *it_branch;
				if(node < N){
					derived_allele[0] = (*it_mut).mutation_type[2];
					ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, 0, derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}else{
					derived_allele[0] = (*it_mut).mutation_type[2];
					assert(prev_rewire[node] == node);
					assert(convert_nodes_prev[node] != 0);
					ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes_prev[node], TSK_NULL, 0, derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}
				site_count++;
			}

			l++;
			it_mut++; 
			if(l == L) break;
		}
		snp_end = l;
		if(snp_end < L){
			pos_end = (*it_mut).pos;
		}else{
			pos_end = (*std::prev(ancmut.mut_end(),1)).pos + 1;
		}

		//write edges

		//node is new if update[n] == 0 and rewire[n] == n
		//if rewire[n] != n, this node is removed
		//if update[n] > 0 and rewire[n] == n, then it is identical to a node in the prev tree
		for(int n = 0; n < 2*N-2; n++){
			int parent_prev = prev_rewire[(*prev_mtr.tree.nodes[n].parent).label];
			int n_now       = update_forwards[n];
			int parent_now  = rewire[(*mtr.tree.nodes[n_now].parent).label]; 
			//rewire might still change, however if it was identical to parent_prev, this would be already reflected here

			if(n < N){
				if( update_forwards[parent_prev] == 0 || update_forwards[parent_prev] != parent_now ){
					//these edges don't exist anymore
					if(n > 0) assert(convert_nodes_prev[n] != 0);
					assert(convert_nodes_prev[parent_prev] != 0);
					ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n], NULL, 0);
					check_tsk_error(ret); 
					edge_count++;
					//prev_mtr.tree.nodes[n].SNP_begin = pos_end;
					mtr.tree.nodes[n].SNP_begin = pos_end; 
				}
			}else if( prev_rewire[n] == n && (n_now == 0 || update_forwards[parent_prev] == 0 || parent_now == 0 || update_forwards[parent_prev] != parent_now) ){
				//these edges don't exist anymore
				assert(convert_nodes_prev[n] != 0);
				assert(convert_nodes_prev[parent_prev] != 0);
				ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n], NULL, 0);
				check_tsk_error(ret);
				if(n_now > 0) mtr.tree.nodes[n_now].SNP_begin = pos_end; 
				edge_count++; 
			}

		}

		//load next tree
		prev_mtr                = mtr;
		prev_leaves             = leaves;
		prev_rewire             = rewire;
		convert_nodes_prev      = convert_nodes;
		it_mut_prev_begin       = it_mut_begin;
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

		//skip last part if I have just read the last tree
		if(num_bases_tree_persists >= 0){    

			it_mut_begin            = it_mut;
			it_mut                  = it_mut_prev_begin;
			std::fill(std::next(rewire.begin(),N), std::prev(rewire.end(),1), 0);
			std::fill(std::next(convert_nodes.begin(),N), std::prev(convert_nodes.end(),1), 0);
			mtr.tree.FindAllLeaves(leaves);

			//remap mutations
			it_mut = it_mut_begin;
			while((*it_mut).tree == ancmut.get_treecount()){
				m_hap.ReadSNP(sequence, bp);
				if(bp != (*it_mut).pos){
					std::cerr << "Error: haps file and anc/mut files don't contain same set of SNPs." << std::endl;
					exit(1);
				}

				sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
				for(int i = 0; i < data.N; i++){
					if(sequence[i] == '1'){
						sequences_carrying_mutation.member[i] = 1;
						sequences_carrying_mutation.num_leaves++;
					}else{
						sequences_carrying_mutation.member[i] = 0;
					}
				}

				if(sequences_carrying_mutation.num_leaves > 0 && sequences_carrying_mutation.num_leaves < N){
					MapMutationExact(mtr.tree, sequences_carrying_mutation, it_mut);
				}
				it_mut++;
			}
			it_mut = it_mut_prev_begin;

			//init
			for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node ++){
				(*it_node).SNP_begin = (*it_mut_begin).pos;
			}
			for(int i = 0; i < 2*N-2; i++){
				if(mtr.tree.nodes[i].num_events >= 1.0){
					rewire[i] = i;
				}
			}

			//////
			//identify equivalent branches,
			//populate update_forwards and update_backwards
			//set prev_rewire[i] = i and rewire[k] = k if these nodes should exist
			//transfer relevant nodes from convert_nodes_prev to convert_nodes 
			FindIdenticalNodes(prev_mtr.tree, mtr.tree, leaves, prev_leaves, update_forwards, update_backwards, prev_rewire, rewire, convert_nodes, convert_nodes_prev, N);
			for(int i = N; i < 2*N-2; i++){
				if(prev_mtr.tree.nodes[i].num_events < 1.0){
					int k = i;
					while(prev_rewire[k] != k){
						k         = (*prev_mtr.tree.nodes[k].parent).label; 
						prev_rewire[i] = k;
						if(k == root){
							prev_rewire[i] = root;
							break;
						}
					}
				}else{
					assert(prev_rewire[i] == i);
				}
			}

			//write new nodes
			for(int n = N; n < 2*N-2; n++){
				if(update_backwards[n] == 0 && rewire[n] == n){
					convert_nodes[n] = node_count;
					//new node, so add to node table
					//need to think how I can replace num_leaves by actual age
					ret = tsk_node_table_add_row(&tables.nodes, 0, leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
					check_tsk_error(ret); 
					mtr.tree.nodes[n].SNP_begin = (*it_mut_begin).pos; 
					node_count++;
				}
				if(prev_rewire[n] == n && convert_nodes_prev[n] == 0){
					convert_nodes_prev[n] = node_count;
					//assert(update_forwards[n] != 0);
					prev_mtr.tree.nodes[n].SNP_begin = (*it_mut_prev_begin).pos;
					assert(update_forwards[n] != 0);
					if(update_forwards[n] != 0){
						convert_nodes[update_forwards[n]] = node_count;
						mtr.tree.nodes[update_forwards[n]].SNP_begin = prev_mtr.tree.nodes[n].SNP_begin;
					}
					ret = tsk_node_table_add_row(&tables.nodes, 0, prev_leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0); 
					check_tsk_error(ret); 
					node_count++;
				}
				if(prev_rewire[n] == n){
					assert(convert_nodes_prev[n] > 0);
				}
				if(convert_nodes_prev[n] != 0){
					assert(prev_rewire[n] == n);
				}
			}      

		}   

	} 


	for(int i = N; i < 2*N-2; i++){
		if(prev_mtr.tree.nodes[i].num_events < 1.0){
			int k = i;
			while(prev_rewire[k] != k && prev_mtr.tree.nodes[k].num_events < 1.0){
				k         = (*prev_mtr.tree.nodes[k].parent).label; 
				prev_rewire[i] = k;
				if(k == root){
					prev_rewire[i] = root;
					break;
				}
			}
		}
	}

	//last tree
	pos = (*it_mut_prev_begin).pos;
	snp = prev_mtr.pos;
	tree_count = (*it_mut).tree;
	//Mutation table
	int l = snp;
	while((*it_mut).tree == tree_count){

		for(std::vector<int>::iterator it_branch = (*it_mut).branch.begin(); it_branch != (*it_mut).branch.end(); it_branch++){ 
			node = *it_branch;
			if(node < N){
				derived_allele[0] = (*it_mut).mutation_type[2];
				ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, TSK_UNKNOWN_TIME, derived_allele, 1, NULL, 0);
				check_tsk_error(ret);
			}else{
				derived_allele[0] = (*it_mut).mutation_type[2];
				assert(prev_rewire[node] == node);
				assert(convert_nodes_prev[node] != 0);
				ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes_prev[node], TSK_NULL, TSK_UNKNOWN_TIME, derived_allele, 1, NULL, 0);
				check_tsk_error(ret);
			}
			site_count++;
		}

		l++;
		it_mut++; 
		if(l == L) break;
	}
	snp_end = l;
	if(snp_end < L){
		pos_end = (*it_mut).pos;
	}else{
		pos_end = (*std::prev(ancmut.mut_end(),1)).pos + 1;
	}

	//for last tree need to dump all edges
	for(int n = 0; n < 2*N-2; n++){
		//Edge table
		//these edges don't exist anymore
		if(rewire[n] == n){

			if(n > 0) assert(convert_nodes_prev[n] != 0);
			int parent_prev = prev_rewire[(*prev_mtr.tree.nodes[n].parent).label];
			assert(convert_nodes_prev[parent_prev] != 0);
			ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n], NULL, 0);   
			check_tsk_error(ret);
			edge_count++;
		}
	}

	std::cerr << "Node count; edge count; tree count" << std::endl;
	std::cerr << node_count << " " << edge_count << " " << tree_count << std::endl;

	ret = tsk_table_collection_sort(&tables, NULL, 0);
	check_tsk_error(ret);
	ret = tsk_table_collection_build_index(&tables, 0);
	check_tsk_error(ret);
	//////////////////////////

	// Write out the tree sequence
	ret = tsk_table_collection_dump(&tables, filename_output.c_str(), 0);        
	check_tsk_error(ret);
	tsk_table_collection_free(&tables); 

}

//convert tree sequence to anc/mut
void
ConvertFromTreeSequence(const std::string& filename_anc, const std::string& filename_mut, const std::string& filename_input, bool no_bl, const int seed = std::time(0) + getpid()){

	int ret, iter;

	std::mt19937 rng;
	rng.seed(seed);
	std::uniform_real_distribution<double> dist(0,1);

	tsk_treeseq_t ts;
	tsk_tree_t tree;
	const tsk_site_t *sites;
	tsk_size_t sites_length;
	tsk_size_t num_SNPs, N, num_trees;
	const tsk_mutation_t *mutation;
	const tsk_id_t *samples;
	tsk_id_t *stack, *node_conversion;
	tsk_id_t u, v, root;
	int j, k, stack_top;

	//load tree sequence
	ret = tsk_treeseq_load(&ts, filename_input.c_str(), 0);
	check_tsk_error(ret);
	//initialise tree
	ret = tsk_tree_init(&tree, &ts, 0);
	check_tsk_error(ret);
	//get sample nodes
	//ret = tsk_treeseq_get_samples(&ts, &samples);
	samples = tsk_treeseq_get_samples(&ts);

	//get num sites
	num_SNPs = tsk_treeseq_get_num_sites(&ts);
	//get num samples
	N = tsk_treeseq_get_num_samples(&ts);
	//get num trees
	num_trees = tsk_treeseq_get_num_trees(&ts);

	stack           = (tsk_id_t *) malloc(tsk_treeseq_get_num_nodes(&ts) * sizeof(*stack));
	if (stack == NULL){
		errx(EXIT_FAILURE, "No memory");
	}
	node_conversion = (tsk_id_t *) malloc(tsk_treeseq_get_num_nodes(&ts) * sizeof(*node_conversion));
	if (node_conversion == NULL){
		errx(EXIT_FAILURE, "No memory");
	}
	for(int i = 0; i < tsk_treeseq_get_num_nodes(&ts); i++){
		node_conversion[i] = -1;
	}

	std::cerr << tsk_treeseq_get_num_individuals(&ts) << " " << N << std::endl;
	j = 0;
	for(int i = 0; i < tsk_treeseq_get_num_individuals(&ts); i++){
		tsk_individual_t *ind;
		ind = (tsk_individual_t *) malloc(2 * sizeof(*ind));
		tsk_treeseq_get_individual(&ts, i, ind);
		//std::cerr << "Node: " << i << std::endl;
		for(int k = 0; k < ind->nodes_length; k++){
			if((ind->nodes[k]) < N){
				node_conversion[j] = (ind->nodes[k]);
				//std::cerr << node_conversion[j] << " ";
				j++;
			}
		}
		//std::cerr << std::endl;
		//std::cerr << ind->nodes[0] << " " << ind->nodes[1] << std::endl;
	}

	//anc/mut variables
	int snp = 0, SNP_begin, SNP_end;
	int bp  = 0, left_bp, right_bp;
	int tree_count = 0;
	int node_count = 0, parent, node, num_children;
	double t1, t2;
	std::string allele;

	MarginalTree mtr;
	Mutations mut;

	mut.info.resize(num_SNPs);

	//count number of trees with at least one SNP
	num_trees = 0;
	//forward iteration through tree sequence
	for(iter = tsk_tree_first(&tree); iter == 1; iter = tsk_tree_next(&tree)){

		//get sites and mutations
		ret = tsk_tree_get_sites(&tree, &sites, &sites_length);
		check_tsk_error(ret);
		//only store tree if it contains at least one site
		bool include = false;
		if(sites_length > 0 && !include){  
			for(j = 0; j < sites_length; j++){
				if(sites[j].mutations_length == 1 && sites[j].ancestral_state_length == 1){
					mutation = &sites[j].mutations[0]; //only one mutation
					if(mutation -> derived_state_length == 1){
						include = true;
						break;
					}
				}
			}
		}
		if(include) num_trees++;

	}

	std::ofstream os(filename_anc);
	FILE *fp = std::fopen(filename_anc.c_str(), "w");
	fprintf(fp, "NUM_HAPLOTYPES %llu ", N);
	tsk_tree_first(&tree);
	std::vector<double> sample_ages(N, 0.0);
	bool any_ancient = false;
	for(int n = 0; n < N; n++){
		tsk_tree_get_time(&tree, n, &sample_ages[node_conversion[n]]);
		if(sample_ages[node_conversion[n]] > 0) any_ancient = true;
	}
	if(any_ancient){
		for(int n = 0; n < N; n++){
			fprintf(fp, "%f ", sample_ages[n]);
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "NUM_TREES %llu\n", num_trees);

	//forward iteration through tree sequence
	for(iter = tsk_tree_first(&tree); iter == 1; iter = tsk_tree_next(&tree)){

		int ntip = 0;
		//get sites and mutations
		ret = tsk_tree_get_sites(&tree, &sites, &sites_length);
		check_tsk_error(ret);

		std::vector<int> used(2*N-1,0);

		//only store tree if it contains at least one site
		bool include = false;
		if(sites_length > 0 && !include){  
			for(j = 0; j < sites_length; j++){
				if(sites[j].mutations_length == 1 && sites[j].ancestral_state_length == 1){
					mutation = &sites[j].mutations[0]; //only one mutation
					if(mutation -> derived_state_length == 1){
						include = true;
						break;
					}
				}
			}
		}
		if(include){

			//for(int i = 0; i < tsk_treeseq_get_num_nodes(&ts); i++){
			//	node_conversion[i] = -1;
			//}

			mtr.pos = snp;
			mtr.tree.nodes.clear();
			mtr.tree.nodes.resize(2*N-1); 

			left_bp = tree.interval.left;
			right_bp = tree.interval.right;

			//get topology of this tree
			if(tsk_tree_get_num_roots(&tree) > 1){
				errx(EXIT_FAILURE, "Multiple roots in tree.");
			}
			root = tsk_tree_get_left_root(&tree);
			//root = tree.left_root;
			node_count = 2*N-2;
			node_conversion[root] = node_count;
			mtr.tree.nodes[node_count].label = node_count;
			used[node_count] = root;
			node_count--;
			stack_top  = 0;
			stack[stack_top] = root;

			//go from root to leaves
			//start with 2N-2, decrease for each subsequent node
			//Need an array saying node x in tree sequence is node y in anc/mut

			while(stack_top >= 0){
				u = stack[stack_top];
				stack_top--;

				//std::cerr << node_count << std::endl;
				num_children = 0;
				for(v = tree.left_child[u]; v != TSK_NULL; v = tree.right_sib[v]) num_children++;

				if(num_children > 0){
					//if(!tsk_treeseq_is_sample(&ts, u)){

					if(tsk_treeseq_is_sample(&ts, u)){
						if(node_conversion[u] == -1){
							assert(u < N);
							node_conversion[u] = u; //TODO
							ntip++;
						}
					}else{
						bool renew = (node_conversion[u] == -1);
						if(node_conversion[u] != -1){
							renew = (renew || (used[node_conversion[u]] != u));
						}
						if(renew){
							node_conversion[u] = node_count;
							used[node_count] = u;
							node_count--;
						}
					}
					parent = node_conversion[u];
					assert(parent != -1);

					if(0){
					if(num_children != 2){
						std::cerr << "nchildren: " << num_children << std::endl;
						for(v = tree.left_child[u]; v != TSK_NULL; v = tree.right_sib[v]){
							if(num_children == 1){
								std::cerr << u << " " << v << " " << "is_sample: " << tsk_treeseq_is_sample(&ts, v) << std::endl;
							}
						}
					}
					}

					if(num_children <= 2){
						for(v = tree.left_child[u]; v != TSK_NULL; v = tree.right_sib[v]){

							if(tsk_treeseq_is_sample(&ts, v)){
								if(node_conversion[v] == -1){
									assert(v < N);
									node_conversion[v] = v; //TODO
									ntip++;
								}
							}else{
								bool renew = (node_conversion[v] == -1);
								if(node_conversion[v] != -1){
									renew = (renew || (used[node_conversion[v]] != v));
								}
								if(renew){
									node_conversion[v] = node_count;
									used[node_count] = v;
									node_count--;
								}
							}
							node = node_conversion[v];

							//if(node == 270) std::cerr << "debug " << node << " " << parent << std::endl;
							mtr.tree.nodes[node].parent    = &mtr.tree.nodes[parent];
							//mtr.tree.nodes[parent].label   = parent;
							mtr.tree.nodes[node].label     = node; 
							if(mtr.tree.nodes[parent].child_left == NULL){
								mtr.tree.nodes[parent].child_left  = &mtr.tree.nodes[node];
							}else{
								mtr.tree.nodes[parent].child_right = &mtr.tree.nodes[node];
							} 

							tsk_tree_get_time(&tree, v, &t1);
							tsk_tree_get_time(&tree, tree.parent[v], &t2);
							mtr.tree.nodes[node].branch_length = t2 - t1;
							mtr.tree.nodes[node].SNP_begin = snp; //first SNP index
						} 
					}else if(num_children > 2){

						//assert(num_children == 2);
						//break polytomies at random
						std::vector<tsk_id_t> children(num_children);
						int i = 0;
						int tmp_node_count, num_children_total = num_children;
						for(v = tree.left_child[u]; v != TSK_NULL; v = tree.right_sib[v]){ 
							if(tsk_treeseq_is_sample(&ts, v)){
								if(node_conversion[v] == -1){
									assert(v < N);
									node_conversion[v] = v; //TODO
									ntip++;
								}
							}else{

								bool renew = (node_conversion[v] == -1);
								if(node_conversion[v] != -1){
									renew = (renew || (used[node_conversion[v]] != v));
								}
								if(renew){
									i++;
								}
							}
						}

						tmp_node_count = node_count - i - (num_children - 2) + 1;
						int min_node_count = tmp_node_count, max_node_count = node_count; //for debugging

						i = 0;
						for(v = tree.left_child[u]; v != TSK_NULL; v = tree.right_sib[v]){ 
							if(!tsk_treeseq_is_sample(&ts, v)){
								bool renew = (node_conversion[v] == -1);
								if(node_conversion[v] != -1){
									renew = (renew || (used[node_conversion[v]] != v));
								}

								if(renew){
									node_conversion[v] = tmp_node_count;
									used[tmp_node_count] = v;
									tmp_node_count++;
									node_count--;
								}
							}
							children[i] = node_conversion[v];
							i++;
						}

						//choose two children at random, in children, replace one of them  
						int childA, childB, ichildA, ichildB, parent_all; //assume randomly chosen
						parent_all = parent;

						tsk_tree_get_time(&tree, tree.left_child[u], &t1);
						tsk_tree_get_time(&tree, u, &t2);

						parent = tmp_node_count;
						tmp_node_count++;
						node_count--;
						while(num_children > 2){
							ichildA = dist(rng)*num_children;
							ichildB = dist(rng)*(num_children - 1);
							if(ichildB >= ichildA) ichildB++;

							childA  = children[ichildA];
							childB  = children[ichildB];

							mtr.tree.nodes[childA].parent    = &mtr.tree.nodes[parent]; 
							mtr.tree.nodes[childA].label     = childA; 
							mtr.tree.nodes[parent].child_left  = &mtr.tree.nodes[childA];
							mtr.tree.nodes[childA].branch_length = 0.0;
							mtr.tree.nodes[childA].SNP_begin = snp; //first SNP index

							mtr.tree.nodes[childB].parent    = &mtr.tree.nodes[parent]; 
							mtr.tree.nodes[childB].label     = childB; 
							mtr.tree.nodes[parent].child_right  = &mtr.tree.nodes[childB];
							mtr.tree.nodes[childB].branch_length = 0.0;
							mtr.tree.nodes[childB].SNP_begin = snp; //first SNP index

							if(ichildA < ichildB){
								children[ichildA] = parent;
								children[ichildB] = children[num_children-1];
								num_children--;
							}else{              
								children[ichildB] = parent;
								children[ichildA] = children[num_children-1];
								num_children--;
							}
							parent = tmp_node_count;
							tmp_node_count++;
							node_count--;
						}
						node_count++;
						tmp_node_count--;
						assert(node_count == min_node_count - 1);
						assert(tmp_node_count == max_node_count + 1);
						parent = parent_all;
						childA = children[0];
						childB = children[1];

						mtr.tree.nodes[childA].parent    = &mtr.tree.nodes[parent]; 
						mtr.tree.nodes[childA].label     = childA; 
						mtr.tree.nodes[parent].child_left  = &mtr.tree.nodes[childA];
						mtr.tree.nodes[childA].branch_length = t2 - t1;
						mtr.tree.nodes[childA].SNP_begin = snp; //first SNP index

						mtr.tree.nodes[childB].parent    = &mtr.tree.nodes[parent]; 
						mtr.tree.nodes[childB].label     = childB; 
						mtr.tree.nodes[parent].child_right  = &mtr.tree.nodes[childB];
						mtr.tree.nodes[childB].branch_length = t2 - t1;
						mtr.tree.nodes[childB].SNP_begin = snp; //first SNP index

					}

				}

				for(v = tree.left_child[u]; v != TSK_NULL; v = tree.right_sib[v]){
					stack_top++;
					stack[stack_top] = v;
				}
				}

				//go through all tip nodes and if they have children, add nodes
				for(int i = 0; i < N; i++){
					//std::cerr << mtr.tree.nodes[i].label << " " << (*mtr.tree.nodes[i].parent).label << std::endl;
					if(mtr.tree.nodes[i].child_left != NULL){
						int child = (*mtr.tree.nodes[i].child_left).label;
						int parent = (*mtr.tree.nodes[i].parent).label;

						mtr.tree.nodes[node_count].label = node_count;
						mtr.tree.nodes[node_count].child_left = &mtr.tree.nodes[i];
						mtr.tree.nodes[node_count].child_right = &mtr.tree.nodes[child];
						mtr.tree.nodes[i].parent = &mtr.tree.nodes[node_count];

						mtr.tree.nodes[node_count].parent = &mtr.tree.nodes[parent];
            if((*mtr.tree.nodes[parent].child_left).label == i){
							mtr.tree.nodes[parent].child_left = &mtr.tree.nodes[node_count];
						}
						if((*mtr.tree.nodes[parent].child_right).label == i){
							mtr.tree.nodes[parent].child_right = &mtr.tree.nodes[node_count];
						}

						mtr.tree.nodes[child].parent = &mtr.tree.nodes[node_count];

						mtr.tree.nodes[node_count].SNP_begin = snp;
						mtr.tree.nodes[node_count].branch_length = mtr.tree.nodes[child].branch_length;
						mtr.tree.nodes[child].branch_length = 0.0;
						mtr.tree.nodes[i].child_left = NULL;
						mtr.tree.nodes[i].child_right = NULL;
						node_count--;
						//std::cerr << i << " " << (*mtr.tree.nodes[i].child_left).label << std::endl;
					}
					if(mtr.tree.nodes[i].child_right != NULL){
						int child = (*mtr.tree.nodes[i].child_right).label;
						int parent = (*mtr.tree.nodes[i].parent).label;

						mtr.tree.nodes[node_count].label = node_count;
						mtr.tree.nodes[node_count].child_left = &mtr.tree.nodes[i];
						mtr.tree.nodes[node_count].child_right = &mtr.tree.nodes[child];
						mtr.tree.nodes[i].parent = &mtr.tree.nodes[node_count];

						mtr.tree.nodes[node_count].parent = &mtr.tree.nodes[parent];
            if((*mtr.tree.nodes[parent].child_left).label == i){
							mtr.tree.nodes[parent].child_left = &mtr.tree.nodes[node_count];
						}
						if((*mtr.tree.nodes[parent].child_right).label == i){
							mtr.tree.nodes[parent].child_right = &mtr.tree.nodes[node_count];
						}

						mtr.tree.nodes[child].parent = &mtr.tree.nodes[node_count];

						mtr.tree.nodes[node_count].SNP_begin = snp;
						mtr.tree.nodes[node_count].branch_length = mtr.tree.nodes[child].branch_length;
						mtr.tree.nodes[child].branch_length = 0.0;
						mtr.tree.nodes[i].child_left = NULL;
						mtr.tree.nodes[i].child_right = NULL;
						node_count--;
						//std::cerr << i << " " << (*mtr.tree.nodes[i].child_right).label << std::endl;
					}
				}
				//std::cerr << ntip << " " << node_count << " " << N-1 << std::endl;

				std::vector<int> checkp(2*N-1,0);
				for(int i = 0; i < 2*N-2; i++){
					assert(mtr.tree.nodes[i].label == i);
          checkp[(*mtr.tree.nodes[i].parent).label]++;
				}
				for(int i = N; i < 2*N-1; i++){
          assert(checkp[i] == 2);
				}

				std::vector<int> checkc(2*N-1,0);
				for(int i = N; i < 2*N-1; i++){
					if(mtr.tree.nodes[i].child_left == NULL) std::cerr << i << std::endl;
					if(mtr.tree.nodes[i].child_right == NULL) std::cerr << i << std::endl;
					checkc[(*mtr.tree.nodes[i].child_left).label]++;
					checkc[(*mtr.tree.nodes[i].child_right).label]++;
				}
				for(int i = 0; i < 2*N-2; i++){
					assert(checkc[i] == 1);
				}
				assert(node_count == N-1); 

				std::vector<float> coords(2*N-1, 0.0);
				if(no_bl){

					int num_lin, parent_num_lin; 
					assert(coords.size() == 2*N-1);
					for(int i = N; i < 2*N-1; i++){
						num_lin = (2*N - i);
						coords[i] = coords[i-1] + 15000.0/(num_lin * (num_lin - 1.0));
					}

					for(int i = 0; i < 2*N-2; i++){
						parent = (*mtr.tree.nodes[i].parent).label;
						mtr.tree.nodes[i].branch_length = (coords[parent] - coords[i]);      
					} 

				}

				for(j = 0; j < sites_length; j++){
					if(sites[j].mutations_length == 1){
						if(sites[j].ancestral_state_length == 1){
							mutation = &sites[j].mutations[0]; //only one mutation
							if(mutation -> derived_state_length == 1){
								//sites[j].pos, mut.id, mut.node, mut.derived_state_length, mut.derived_state
								mut.info[snp].snp_id = snp;
								mut.info[snp].pos    = round(sites[j].position);
								mut.info[snp].tree   = tree_count;

								allele.assign(sites[j].ancestral_state, sites[j].ancestral_state_length);
								mut.info[snp].mutation_type = allele + "/";
								allele.assign(mutation -> derived_state, mutation -> derived_state_length);
								mut.info[snp].mutation_type += allele;

								mut.info[snp].branch.resize(1);
								node                    = node_conversion[mutation -> node];
								mut.info[snp].branch[0] = node;
								mut.info[snp].rs_id     = std::to_string(mutation -> id);
								if(!no_bl){
									tsk_tree_get_time(&tree, mutation -> node, &t1);
									tsk_tree_get_time(&tree, tree.parent[mutation -> node], &t2);
									mut.info[snp].age_begin = t1;
									mut.info[snp].age_end   = t2;
								}else{
									mut.info[snp].age_begin = coords[node];
									mut.info[snp].age_end   = coords[(*mtr.tree.nodes[node].parent).label];             
								}
								mtr.tree.nodes[node].num_events += 1.0;

								if(snp > 0){
									mut.info[snp-1].dist = mut.info[snp].pos - mut.info[snp-1].pos;
								}
								snp++;
							}
						}
					}
				}

				SNP_end = snp-1;
				mut.info[SNP_end].dist = 1.0;
				for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
					(*it_node).SNP_end = SNP_end; //last SNP index 
				}

				mtr.Dump(fp);
				tree_count++; 
			}

		}

		std::fclose(fp);

		//Dump mut file
		mut.info.resize(snp);
		mut.header = "snp;pos_of_snp;dist;rs-id;tree_index;branch_indices;is_not_mapping;is_flipped;age_begin;age_end;ancestral_allele/alternative_allele;";
		mut.Dump(filename_mut);

	}


#endif //TREE_SEQUENCE_HPP 
