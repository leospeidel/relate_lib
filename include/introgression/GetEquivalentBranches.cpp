#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "anc.hpp"
#include "mutations.hpp"
#include "cxxopts.hpp"
#include <ctime>

struct branch{
  Leaves leaves, leaves_parent;
  SNPInfo snp_info; 
  int SNP_begin, SNP_end; 
  std::vector<double> lower_age, upper_age;
  int num_events;
};

struct equiv_branches{
  SNPInfo snp_info; 
  double lower_age, upper_age;
  int num_events;
  std::vector<int> bp;
  std::vector<float> w1;
};

void
GetEquivalentBranches(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Get list of mutations on branches." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Get list of mutations on branches ..." << std::endl;


  //For branches of interest, determine how many of them are equivalent.
  //For equivalent branches, calculate mean of lower age and upper age
  //Output as .mut file

  std::string line;
  //////////////////////////////////
  //Parse Data
  int N;
  igzstream is_N(options["input"].as<std::string>() + ".anc.gz");
  if(is_N.fail()){
    std::cerr << "Error while reading anc." << std::endl;
    exit(1);
  }
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  int N_total = 2*N-1;
  Data data(N, 1);

  std::cerr << N << std::endl;
 

  //////////////////////////////////////////// Read Tree ///////////////////////////////////

  //read anc
  AncesTree anc;
  anc.Read(options["input"].as<std::string>() + ".anc");
  //read mut
  Mutations mut;
  mut.Read(options["input"].as<std::string>() + ".mut");

  //Associate branches
  //Pre calculate how many descendants a branch needs to be equivalent
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

  //find equivalent branches
  CorrTrees::iterator it_seq_prev;
  CorrTrees::iterator it_seq; 
  it_seq_prev = anc.seq.begin();
  it_seq      = std::next(it_seq_prev,1); 

  std::vector<std::vector<int>> equivalent_branches;
  std::vector<std::vector<int>>::iterator it_equivalent_branches;
  std::vector<std::vector<int>>::reverse_iterator rit_equivalent_branches;

  for(; it_seq != anc.seq.end();){
    equivalent_branches.emplace_back();
    it_equivalent_branches = std::prev(equivalent_branches.end(),1);
    anc.BranchAssociation((*it_seq_prev).tree, (*it_seq).tree, *it_equivalent_branches, potential_branches, N, N_total, threshold_brancheq); //O(N^2) 
    it_seq++;
    it_seq_prev++;
  }  
 
  equivalent_branches.emplace_back();
  it_equivalent_branches = std::prev(equivalent_branches.end(),1);
  (*it_equivalent_branches).resize(2*N-1);
  for(std::vector<int>::iterator it_b = (*it_equivalent_branches).begin(); it_b != (*it_equivalent_branches).end(); it_b++){
    *it_b = -1;
  }

  int num_trees = equivalent_branches.size();

  //use equivalent branches as follows:
  //read in mutation, treeID, branchID. follow equivalent branches to get all mutations on this branch
  //mark these mutations as read.
  std::vector<equiv_branches> ebranches;
  Mutations ebranches_mut;
  
  int L = mut.info.size();
  std::vector<int> checked(L, 0);
 
  for(int snp = mut.info.size()-1; snp >= 1; snp--){
 
    if(checked[snp] == 0 && mut.info[snp].branch.size() == 1){

      checked[snp] = 1;

      equiv_branches b;
      ebranches_mut.info.push_back(mut.info[snp]);
      b.bp.push_back(mut.info[snp].pos);

      int treeID   = mut.info[snp].tree;
      int branchID = mut.info[snp].branch[0];
      int snp_tmp = snp-1;
      while(mut.info[snp_tmp].tree == treeID && snp_tmp >= 0){
        if(mut.info[snp_tmp].branch.size() == 1 && checked[snp_tmp] == 0){
          if(mut.info[snp_tmp].branch[0] == branchID){
            b.bp.push_back(mut.info[snp_tmp].pos);
            b.lower_age += mut.info[snp_tmp].age_begin;
            b.upper_age += mut.info[snp_tmp].age_end;
            checked[snp_tmp] = 1;
          }
        }
        snp_tmp--;
        if(snp_tmp == -1) break;
      }

      treeID--;
      if(treeID >= 0){
        branchID = equivalent_branches[treeID][branchID];
      }else{
        branchID = -1;
      }

      while(branchID != -1 && treeID >= 0){
       
        while(mut.info[snp_tmp].tree == treeID && snp_tmp >= 0){
          if(mut.info[snp_tmp].branch.size() == 1 && checked[snp_tmp] == 0){
            if(mut.info[snp_tmp].branch[0] == branchID){
              b.bp.push_back(mut.info[snp_tmp].pos);
              b.lower_age += mut.info[snp_tmp].age_begin;
              b.upper_age += mut.info[snp_tmp].age_end;
              checked[snp_tmp] = 1;
            }
          }
          snp_tmp--;
          if(snp_tmp == -1) break;
        }

        if(snp_tmp == -1) break;
        treeID--;
        if(treeID >= 0){
          branchID = equivalent_branches[treeID][branchID];
        }else{
          branchID = -1;
        }

      }

      ebranches.push_back(b);
    
    }
   
  }

  std::ofstream os(options["output"].as<std::string>() + ".equiv.mut.map");

  os << "treeID branchID bp\n";

  for(int snp = 0; snp < ebranches.size(); snp++){
  
    for(std::vector<int>::iterator it_bp = ebranches[snp].bp.begin(); it_bp != ebranches[snp].bp.end(); it_bp++){
      os << ebranches_mut.info[snp].tree << " " << ebranches_mut.info[snp].branch[0] << " " << *it_bp << "\n";
    }
  
  }
  os.close();

  ebranches_mut.Dump(options["output"].as<std::string>() + ".equiv.mut");

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}

