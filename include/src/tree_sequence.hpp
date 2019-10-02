#ifndef TREE_SEQUENCE_HPP
#define TREE_SEQUENCE_HPP

#include <iostream>
#include <err.h>
#include <deque>

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

void
PropagateMutationExact(Node& node, std::deque<int>& branches, std::deque<int>& branches_flipped, Leaves& sequences_carrying_mutations, PropagateStructLocal& report){

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
  std::deque<int> branches;
  std::deque<int> branches_flipped; 

  PropagateStructLocal report;
  PropagateMutationExact(*std::prev(tree.nodes.end(), 1), branches, branches_flipped, sequences_carrying_mutations, report);

  if( branches_flipped.size() == 0 ){

    assert(branches.size() > 0);
    (*it_mut).branch  = branches;
    for(std::deque<int>::iterator it = branches.begin(); it != branches.end(); it++){
      tree.nodes[*it].num_events += 1.0;
    }
    return branches.size();

  }else{

    if( branches.size() <= branches_flipped.size() && branches.size() > 0 ){ 

      (*it_mut).branch  = branches;
      for(std::deque<int>::iterator it = branches.begin(); it != branches.end(); it++){
        tree.nodes[*it].num_events += 1.0;
      }
      return branches.size();

    }else{

      (*it_mut).flipped = true;
      (*it_mut).branch  = branches_flipped;
      for(std::deque<int>::iterator it = branches_flipped.begin(); it != branches_flipped.end(); it++){
        tree.nodes[*it].num_events += 1.0;
      }
      return branches_flipped.size();

    }

  }

}

/////////////////////////////////

//convert to tree sequence (new set of nodes for each tree)
void
DumpAsTreeSequence(const std::string& filename_anc, const std::string& filename_mut, const std::string& filename_output){

  ////////////////////////
  //read in anc file

  MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;

  ////////// 1. Read one tree at a time /////////

  //We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
  //The mut file is read once, file is closed after constructor is called.
  AncMutIterators ancmut(filename_anc, filename_mut);


  num_bases_tree_persists = ancmut.FirstSNP(mtr, it_mut);
  int N = (mtr.tree.nodes.size() + 1)/2.0, root = 2*N - 2, L = ancmut.NumSnps();
  Data data(N,L);
  std::vector<float> coordinates(2*data.N-1,0.0);
  mtr.tree.GetCoordinates(coordinates);

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
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0 , NULL, 0);
  }

  //population table

  //sites table
  char ancestral_allele[1];
  //tsk_site_table_add_row(&tables.sites, 1, ancestral_allele, sizeof(ancestral_allele), NULL, 0);
  for(; it_mut != ancmut.mut_end(); it_mut++){
    ancestral_allele[0] = (*it_mut).mutation_type[0];
    ret = tsk_site_table_add_row(&tables.sites, (*it_mut).pos, ancestral_allele, 1, NULL, 0);
    check_tsk_error(ret);
  }

  for(int i = 0; i < data.N; i++){
    ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, i, NULL, 0);
    check_tsk_error(ret);
  }

  ///////////////////////////////////////////////////////////////////////////// 

  int pos, snp, pos_end, snp_end, tree_count = 0, node, node_const, site_count = 0;
  int node_count = data.N, edge_count = 0;

  char derived_allele[1];
  while(num_bases_tree_persists >= 0.0){

    pos = (*it_mut).pos;
    if(mtr.pos == 0) pos = 0;
    snp = mtr.pos;

    tree_count = (*it_mut).tree;
    node_const = tree_count * (data.N - 1);

    //Mutation table
    int l = snp;
    while((*it_mut).tree == tree_count){
      if((*it_mut).branch.size() == 1){
        node = *(*it_mut).branch.begin();
        if(node < N){
          derived_allele[0] = (*it_mut).mutation_type[2];
          ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, derived_allele, 1, NULL, 0);
          check_tsk_error(ret);
        }else{
          derived_allele[0] = (*it_mut).mutation_type[2];
          ret = tsk_mutation_table_add_row(&tables.mutations, l, node + node_const, TSK_NULL, derived_allele, 1, NULL, 0);
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

    //Node table
    std::vector<Node>::iterator it_node = std::next(mtr.tree.nodes.begin(), data.N);
    int n = N;
    for(std::vector<float>::iterator it_coords = std::next(coordinates.begin(), data.N); it_coords != coordinates.end(); it_coords++){   
      ret = tsk_node_table_add_row(&tables.nodes, 0, *it_coords, TSK_NULL, TSK_NULL, NULL, 0);   
      check_tsk_error(ret);
      n++;
      node_count++;
    }

    //Edge table
    for(it_node = mtr.tree.nodes.begin(); it_node != std::prev(mtr.tree.nodes.end(),1); it_node++){
      node = (*it_node).label;
      if(node >= data.N) node += node_const;
      ret = tsk_edge_table_add_row(&tables.edges, pos, pos_end, (*(*it_node).parent).label + node_const, node);    
      check_tsk_error(ret);
      edge_count++;
    }

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    mtr.tree.GetCoordinates(coordinates);

  } 

  tsk_table_collection_sort(&tables, NULL, 0);
  check_tsk_error(ret);

  std::cerr << "Node count; edge count; tree count" << std::endl;
  std::cerr << node_count << " " << edge_count << " " << tree_count << std::endl;


  //////////////////////////

  // Write out the tree sequence
  ret = tsk_table_collection_dump(&tables, filename_output.c_str(), 0);        
  check_tsk_error(ret);
  tsk_table_collection_free(&tables); 

}

//compress by combining equivalent branches (ignores branch lengths)
void
DumpAsTreeSequenceTopoOnly(const std::string& filename_anc, const std::string& filename_mut, const std::string& filename_output){

  MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  std::vector<Leaves> leaves, prev_leaves;
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;

  ////////// 1. Read one tree at a time /////////

  //We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
  //The mut file is read once, file is closed after constructor is called.
  AncMutIterators ancmut(filename_anc, filename_mut);

  num_bases_tree_persists = ancmut.FirstSNP(mtr, it_mut);
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
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0 , NULL, 0);
  }

  //population table

  //sites table
  char ancestral_allele[1];
  //tsk_site_table_add_row(&tables.sites, 1, ancestral_allele, sizeof(ancestral_allele), NULL, 0);
  for(; it_mut != ancmut.mut_end(); it_mut++){
    ancestral_allele[0] = (*it_mut).mutation_type[0];
    ret = tsk_site_table_add_row(&tables.sites, (*it_mut).pos, ancestral_allele, 1, NULL, 0);
    check_tsk_error(ret);
  }

  //........................................................................
  //Iterate through ancmut

  it_mut = ancmut.mut_begin();
  //std::vector<float> coordinates(2*N-1,0.0);
  int pos, snp, pos_end, snp_end, tree_count = 0, node, site_count = 0;

  int node_count = 0, edge_count = 0;
  bool is_different = false;
  //for each tree, keep a vector convert_nodes that maps nodes to ts nodes
  std::vector<int> convert_nodes(mtr.tree.nodes.size(), 0), convert_nodes_prev(mtr.tree.nodes.size(), 0);
  std::vector<int> update_backwards(2*N-1,0), update_forwards(2*N-1,0);

  std::vector<int>::iterator it_update_backwards = update_backwards.begin(), it_update_forwards = update_forwards.begin();
  std::vector<int>::iterator it_convert = convert_nodes.begin();
  for(; it_convert != std::next(convert_nodes.begin(),N); it_convert++){
    *it_convert = node_count;
    *it_update_backwards = node_count;
    *it_update_forwards  = node_count;
    //new node, so add to node table
    //need to think how I can replace num_leaves by actual age
    ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, leaves[node_count].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);   
    node_count++;
    it_update_forwards++;
    it_update_backwards++;
  }
  for(;it_convert != convert_nodes.end(); it_convert++){
    *it_convert          = node_count;
    //new node, so add to node table
    //need to think how I can replace num_leaves by actual age
    ret = tsk_node_table_add_row(&tables.nodes, 0, leaves[node_count].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);   
    node_count++;
  }

  char derived_allele[1];
  while(num_bases_tree_persists >= 0.0){

    //mtr.tree.GetCoordinates(coordinates);
    pos = (*it_mut).pos;
    if(mtr.pos == 0) pos = 0;
    for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node ++){
      (*it_node).SNP_begin = pos;
    }
    snp = mtr.pos;

    tree_count = (*it_mut).tree;

    if(tree_count > 0){
      //nodes:
      //for each node, check if its descendant set is identical to before
      //if no, update convert_nodes[i] = node_count;
      std::fill(std::next(update_backwards.begin(),N), update_backwards.end(), 0);
      std::fill(std::next(update_forwards.begin(),N), update_forwards.end(), 0);
      std::fill(std::next(convert_nodes_prev.begin(),N), convert_nodes_prev.end(), 0);
      for(int n = 0; n < N; n++){
        mtr.tree.nodes[n].SNP_begin = prev_mtr.tree.nodes[n].SNP_begin;
      }
      //identify all nodes that are new
      for(int n = N; n < 2*N-1; n++){
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

      for(int n = 0; n < 2*N-2; n++){
        int parent_prev = (*prev_mtr.tree.nodes[n].parent).label;
        int n_now       = update_forwards[n];
        int parent_now  = (*mtr.tree.nodes[n_now].parent).label;

        if(n < N){
          if( update_forwards[parent_prev] != parent_now ){
            //these edges don't exist anymore 
            ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes[parent_prev], convert_nodes[n]);
            check_tsk_error(ret); 
            edge_count++;
            mtr.tree.nodes[n].SNP_begin = pos_end; 
          }
        }else if( n_now == 0 || update_forwards[parent_prev] != parent_now ){
          //these edges don't exist anymore
          ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes[parent_prev], convert_nodes[n]);
          if(n_now > 0) mtr.tree.nodes[n_now].SNP_begin = pos_end; 
          check_tsk_error(ret); 
          edge_count++; 
        }
      }

      for(int n = N; n < 2*N-1; n++){
        if(update_backwards[n] == 0){
          convert_nodes[n] = node_count;
          //new node, so add to node table
          //need to think how I can replace num_leaves by actual age
          ret = tsk_node_table_add_row(&tables.nodes, 0, leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);  
          mtr.tree.nodes[n].SNP_begin = pos; 
          node_count++;
        }else{
          convert_nodes[n] = convert_nodes_prev[n];
        }
      }          

    }

    //Mutation table
    int l = snp;
    while((*it_mut).tree == tree_count){
      if((*it_mut).branch.size() == 1){
        node = *(*it_mut).branch.begin();
        if(node < N){
          derived_allele[0] = (*it_mut).mutation_type[2];
          ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, derived_allele, 1, NULL, 0);
          check_tsk_error(ret);
        }else{
          derived_allele[0] = (*it_mut).mutation_type[2];
          ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes[node], TSK_NULL, derived_allele, 1, NULL, 0);
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

    //load next tree
    prev_mtr                = mtr;
    prev_leaves             = leaves;
    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    mtr.tree.FindAllLeaves(leaves);
  } 

  //for last tree need to dump all edges
  for(int n = 0; n < 2*N-2; n++){        
    int parent_prev = (*prev_mtr.tree.nodes[n].parent).label;
    ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes[parent_prev], convert_nodes[n]);
    check_tsk_error(ret); 
    edge_count++;
  }

  std::cerr << "Node count; edge count; tree count" << std::endl;
  std::cerr << node_count << " " << edge_count << " " << tree_count << std::endl;

  tsk_table_collection_sort(&tables, NULL, 0);
  check_tsk_error(ret);

  //////////////////////////

  // Write out the tree sequence
  ret = tsk_table_collection_dump(&tables, filename_output.c_str(), 0);        
  check_tsk_error(ret);
  tsk_table_collection_free(&tables); 

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

//removes branches with no mutation mapped to it (TODO: this is not yet optimal)
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
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0 , NULL, 0);
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
          ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, derived_allele, 1, NULL, 0);
          check_tsk_error(ret);
        }else{
          derived_allele[0] = (*it_mut).mutation_type[2];
          assert(prev_rewire[node] == node);
          assert(convert_nodes_prev[node] != 0);
          ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes_prev[node], TSK_NULL, derived_allele, 1, NULL, 0);
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
          ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n]);
          check_tsk_error(ret); 
          edge_count++;
          //prev_mtr.tree.nodes[n].SNP_begin = pos_end;
          mtr.tree.nodes[n].SNP_begin = pos_end; 
        }
      }else if( prev_rewire[n] == n && (n_now == 0 || update_forwards[parent_prev] == 0 || parent_now == 0 || update_forwards[parent_prev] != parent_now) ){
        //these edges don't exist anymore
        assert(convert_nodes_prev[n] != 0);
        assert(convert_nodes_prev[parent_prev] != 0);
        ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n]);
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
        ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, derived_allele, 1, NULL, 0);
        check_tsk_error(ret);
      }else{
        derived_allele[0] = (*it_mut).mutation_type[2];
        assert(prev_rewire[node] == node);
        assert(convert_nodes_prev[node] != 0);
        ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes_prev[node], TSK_NULL, derived_allele, 1, NULL, 0);
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
      ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n]);   
      check_tsk_error(ret);
      edge_count++;
    }
  }

  std::cerr << "Node count; edge count; tree count" << std::endl;
  std::cerr << node_count << " " << edge_count << " " << tree_count << std::endl;

  tsk_table_collection_sort(&tables, NULL, 0);
  check_tsk_error(ret);

  //////////////////////////

  // Write out the tree sequence
  ret = tsk_table_collection_dump(&tables, filename_output.c_str(), 0);        
  check_tsk_error(ret);
  tsk_table_collection_free(&tables); 

}

//removes branches with no mutation mapped to it, where mutations are remapped so that data can be recovered exactly (TODO: this is not yet optimal)
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
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0 , NULL, 0);
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
      for(std::deque<int>::iterator it_branch = (*it_mut).branch.begin(); it_branch != (*it_mut).branch.end(); it_branch++){ 
        node = *it_branch;
        if(node < N){
          derived_allele[0] = (*it_mut).mutation_type[2];
          ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, derived_allele, 1, NULL, 0);
          check_tsk_error(ret);
        }else{
          derived_allele[0] = (*it_mut).mutation_type[2];
          assert(prev_rewire[node] == node);
          assert(convert_nodes_prev[node] != 0);
          ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes_prev[node], TSK_NULL, derived_allele, 1, NULL, 0);
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
          ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n]);
          check_tsk_error(ret); 
          edge_count++;
          //prev_mtr.tree.nodes[n].SNP_begin = pos_end;
          mtr.tree.nodes[n].SNP_begin = pos_end; 
        }
      }else if( prev_rewire[n] == n && (n_now == 0 || update_forwards[parent_prev] == 0 || parent_now == 0 || update_forwards[parent_prev] != parent_now) ){
        //these edges don't exist anymore
        assert(convert_nodes_prev[n] != 0);
        assert(convert_nodes_prev[parent_prev] != 0);
        ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n]);
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

    for(std::deque<int>::iterator it_branch = (*it_mut).branch.begin(); it_branch != (*it_mut).branch.end(); it_branch++){ 
      node = *it_branch;
      if(node < N){
        derived_allele[0] = (*it_mut).mutation_type[2];
        ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, derived_allele, 1, NULL, 0);
        check_tsk_error(ret);
      }else{
        derived_allele[0] = (*it_mut).mutation_type[2];
        assert(prev_rewire[node] == node);
        assert(convert_nodes_prev[node] != 0);
        ret = tsk_mutation_table_add_row(&tables.mutations, l, convert_nodes_prev[node], TSK_NULL, derived_allele, 1, NULL, 0);
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
      ret = tsk_edge_table_add_row(&tables.edges, prev_mtr.tree.nodes[n].SNP_begin, pos_end, convert_nodes_prev[parent_prev], convert_nodes_prev[n]);   
      check_tsk_error(ret);
      edge_count++;
    }
  }

  std::cerr << "Node count; edge count; tree count" << std::endl;
  std::cerr << node_count << " " << edge_count << " " << tree_count << std::endl;

  tsk_table_collection_sort(&tables, NULL, 0);
  check_tsk_error(ret);

  //////////////////////////

  // Write out the tree sequence
  ret = tsk_table_collection_dump(&tables, filename_output.c_str(), 0);        
  check_tsk_error(ret);
  tsk_table_collection_free(&tables); 

}


#endif //TREE_SEQUENCE_HPP 
