#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "tree_sequence.hpp"
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

void
TMRCAbyTopo(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("haps") || !options.count("sample") || !options.count("input") || !options.count("output") || !options.count("poplabels")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output, poplabels." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate TMRCA by topology of pop_of_interest." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "TMRCA between groups ..." << std::endl;

  std::vector<std::string> filenames_anc, filenames_mut, filenames_chr, filenames_mask, filenames_haps, filenames_sample;

  if(options.count("chr")){
    igzstream is_chr(options["chr"].as<std::string>());
    if(is_chr.fail()){
      std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
    }
    std::string line;
    while(getline(is_chr, line)){
      filenames_anc.push_back(options["input"].as<std::string>() + "_chr" + line + ".anc"); 
      filenames_mut.push_back(options["input"].as<std::string>() + "_chr" + line + ".mut");    
      filenames_haps.push_back(options["haps"].as<std::string>() + "_chr" + line + ".haps.gz"); 
      filenames_sample.push_back(options["sample"].as<std::string>() + "_chr" + line + ".sample.gz"); 
      filenames_chr.push_back(line);
      if(options.count("mask")) filenames_mask.push_back(options["mask"].as<std::string>() + "_chr" + line + ".fa.gz");
    }
    is_chr.close();
  }else{
    filenames_anc.push_back(options["input"].as<std::string>() + ".anc"); 
    filenames_mut.push_back(options["input"].as<std::string>() + ".mut"); 
    filenames_haps.push_back(options["haps"].as<std::string>()); 
    filenames_sample.push_back(options["sample"].as<std::string>()); 
    filenames_chr.push_back("NA");  
    if(options.count("mask")) filenames_mask.push_back(options["mask"].as<std::string>());
  }

  Sample sample;
  sample.Read(options["poplabels"].as<std::string>());
  std::string label;

  if(!options.count("pop_of_interest")){
    label = sample.AssignPopOfInterest("All");
  }else{
    label = sample.AssignPopOfInterest(options["pop_of_interest"].as<std::string>());
  }
  assert(sample.group_of_interest.size() == 3);

  AncMutIterators ancmut_tmp(filenames_anc[0], filenames_mut[0]); 
  Data data(ancmut_tmp.NumTips(),ancmut_tmp.NumSnps());
  assert(data.N % 2 == 0);

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;

  std::vector<int> index;
  //choosing the first haplotype associated with group of interest
  for(int k = 0; k < sample.group_of_interest.size(); k++){
    for(int i = 0; i < data.N; i++){
      if(sample.group_of_haplotype[i] == sample.group_of_interest[k]){
        index.push_back(i);
        break;
      }
    }
  }
  std::sort(index.begin(), index.end());

  std::vector<std::vector<double>> pairwise_tmrca(index.size());
  for(std::vector<std::vector<double>>::iterator it_pairwise_tmrca = pairwise_tmrca.begin(); it_pairwise_tmrca != pairwise_tmrca.end(); it_pairwise_tmrca++){
    (*it_pairwise_tmrca).resize(index.size());
    std::fill((*it_pairwise_tmrca).begin(), (*it_pairwise_tmrca).end(), 0.0);
  }

  ////////// 1. Read one tree at a time /////////

  //We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
  //The mut file is read once, file is closed after constructor is called.

  double genome_length = 0.0;
  std::vector<float> coordinates;
  std::vector<Leaves> desc;

  std::ofstream os_loc_tmrca(options["output"].as<std::string>() + ".topotmrca");
  if(os_loc_tmrca.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".topotmrca" << "." << std::endl;
    exit(1);
  }

  std::ofstream os_mut(options["output"].as<std::string>() + ".topomuts");
  if(os_mut.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".topomuts" << "." << std::endl;
    exit(1);
  }


  os_loc_tmrca << "CHR BP BP_persist num_snps_on_tree frac_branches_with_snp topology";
  for(int i = 0; i < index.size()-1; i++){
    os_loc_tmrca << " tmrca" << i+1;
  }
  os_loc_tmrca << "\n";

  std::vector<int> mut_left, mut_right, mut_both, mut_out;
  int count = 0;
  for(int chr = 0; chr < filenames_anc.size(); chr++){

    std::cerr << "chr: " << chr << std::endl;

    bool first_tree = true;

    //open genealogies
    AncMutIterators ancmut(filenames_anc[chr], filenames_mut[chr]); 
    assert(data.N == ancmut.NumTips());

    //open haps sample files
    haps m_hap(filenames_haps[chr].c_str(), filenames_sample[chr].c_str());
    if(data.N != ancmut.NumTips()){
      std::cerr << "Haps file and anc/mut have different number of samples" << std::endl;
      exit(1);
    }
    int num_snps = ancmut.NumSnps();

    /////////////////////////////////
    //Read anc and mut

    //read in the first 'len' SNPs in haps/sample
    int len = 1000000;
    std::vector<std::vector<char>> sequence(len);
    for(std::vector<std::vector<char>>::iterator it_seq = sequence.begin(); it_seq != sequence.end(); it_seq++){
      (*it_seq).resize(data.N);
    }
    std::vector<int> bp(len,0);
    int snp_mod = 0, snp_index = 0;
    for(snp_mod = 0; snp_mod < std::min(len, m_hap.GetL()); snp_mod++){
      m_hap.ReadSNP(sequence[snp_mod], bp[snp_mod]);
      snp_mod++;
    }
    snp_mod--;

    Leaves sequences_carrying_mutation;
    sequences_carrying_mutation.member.resize(data.N);

    /////////////////////////////////
    //Read anc and mut

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    while(num_bases_tree_persists >= 0.0){

      //new tree

      //calculate some tree properties that correlate with tree quality
      float frac_branches_with_snp = 0, num_snps_on_tree = 0;
      for(std::vector<Node>::iterator it_n = std::next(mtr.tree.nodes.begin(), data.N); it_n != mtr.tree.nodes.end(); it_n++){
        num_snps_on_tree += (*it_n).num_events;
        if( ((*it_n).num_events) >= 1.0) frac_branches_with_snp += 1.0;
      }
      frac_branches_with_snp /= data.N-1.0;

      //by default, I know how many mutations map to each branch in the tree, but I'm not recording which specific mutations these are
      //So I need to remap these.

      //find the maximum range of the tree, i.e. each branch has a range in the genome, find maximum span
      int min_snp = mtr.tree.nodes[0].SNP_begin, max_snp = mtr.tree.nodes[0].SNP_end;
      for(int i = data.N; i < mtr.tree.nodes.size(); i++){
        if(min_snp > mtr.tree.nodes[i].SNP_begin) min_snp = mtr.tree.nodes[i].SNP_begin;
        if(max_snp < mtr.tree.nodes[i].SNP_end) max_snp = mtr.tree.nodes[i].SNP_end;
      }
      int min_bp = ancmut.mut.info[min_snp].pos, max_bp = ancmut.mut.info[max_snp].pos;

      //read all mutations between min_bp and max_bp
      int snp_start = -1, snp_end = -1, bp_dummy = bp[0];
      //bp[snp_mod] is last read entry, i.e. largest BP

      if(max_bp > bp[snp_mod]){ //need to parse additional sites in haps/sample
        while(max_bp > bp[snp_mod]){
          snp_mod++;
          snp_mod = snp_mod % len;
          m_hap.ReadSNP(sequence[snp_mod], bp[snp_mod]);
        }
      }

      //check if first entry is already exceeding
      if(bp[ (snp_mod + 1) % len ] >= min_bp){
        snp_start = (snp_mod + 1) % len;
        if(!first_tree) std::cerr << "Warnings: Haps/sample buffer may be too short." << std::endl;
      }else{
        snp_start = (snp_mod + 1) % len;
        for(int i = 0; i < len; i++){
          if(bp[ (i + snp_mod + 1) % len ] >= min_bp){
            snp_start = (i + snp_mod + 1) % len;
            break;
          }
        }
      }

      if(bp[ (snp_mod + 1) % len ] >= max_bp){
        snp_end = (snp_mod + 1) % len;
        std::cerr << "Warnings: Haps/sample buffer may be too short." << std::endl;
      }else{
        snp_end = (snp_mod + 1) % len;
        for(int i = 0; i < len; i++){
          if(bp[ (i + snp_mod + 1) % len ] >= max_bp){
            snp_end = (i + snp_mod + 1) % len;
            break;
          }
        }
      }

      //need to check for mutations between snp_start and snp_end
      //for each check whether they are mapping
      Tree tr = mtr.tree;
      std::vector<std::vector<int>> mut_on_branches(mtr.tree.nodes.size());
      int num_carriers;	

      //remap mutations
      //get snp index from (*it_mut)
      snp_index = 0;
      //map_mutations
      for(int snp = snp_start; snp <= snp_end; snp++){

        while(ancmut.mut.info[snp_index].pos < bp[snp]) snp_index++;
        if(ancmut.mut.info[snp_index].pos == bp[snp]){

          sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
          for(int i = 0; i < data.N; i++){
            if(sequence[snp % len][i] == '1'){
              sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
              sequences_carrying_mutation.num_leaves++;
            }else{
              sequences_carrying_mutation.member[i] = 0;
            }
          }

          if(sequences_carrying_mutation.num_leaves > 1){
            if(sequences_carrying_mutation.num_leaves > 0 && sequences_carrying_mutation.num_leaves < data.N){

              std::vector<SNPInfo> muts(1);
              Muts::iterator it_m = muts.begin();
              if(MapMutation(tr, sequences_carrying_mutation, it_m, false) < 3){
                int branch = (*it_m).branch[0];
                if(ancmut.mut.info[mtr.tree.nodes[branch].SNP_begin].pos <= bp[snp % len] && ancmut.mut.info[mtr.tree.nodes[branch].SNP_end].pos >= bp[snp % len]){
                  mut_on_branches[branch].push_back(snp_index);
                }
              }

            }
          }

          snp_index++;
        }

      }

      //remapping done


      //now calculate TMRCAs in the tree
      for(std::vector<std::vector<double>>::iterator it_pairwise_tmrca = pairwise_tmrca.begin(); it_pairwise_tmrca != pairwise_tmrca.end(); it_pairwise_tmrca++){
        std::fill((*it_pairwise_tmrca).begin(), (*it_pairwise_tmrca).end(), 0.0);
      }

      int tree_index = (*it_mut).tree;
      genome_length += num_bases_tree_persists/1e12;

      mtr.tree.GetCoordinates(coordinates);
      mtr.tree.FindAllLeaves(desc);

      //check topology
      //find tmrca of all three

      for(int i = data.N; i < 2*data.N-1; i++){

        int child_left  = (*mtr.tree.nodes[i].child_left).label;
        int child_right = (*mtr.tree.nodes[i].child_right).label;
        float coord     = coordinates[i];

        bool any_left = false, any_right = false; //are there any samples in index to the left or right of this coalescence?

        std::vector<bool> is_left(index.size(), false);
        std::vector<bool> is_right(index.size(), false);
        for(int k = 0; k < index.size(); k++){
          for(std::vector<int>::iterator it_hap1 = desc[child_left].member.begin(); it_hap1 != desc[child_left].member.end(); it_hap1++){
            if(*it_hap1 == index[k]){
              is_left[k] = true;
              any_left = true;
            }
          }
        }

        if(any_left){
          for(int k = 0; k < index.size(); k++){
            for(std::vector<int>::iterator it_hap1 = desc[child_right].member.begin(); it_hap1 != desc[child_right].member.end(); it_hap1++){
              if(*it_hap1 == index[k]){
                is_right[k] = true;
                any_right = true;
              }
            }
          }
        }

        //if an individual in index is to either side, this coal event is their MRCA
        if(any_left && any_right){
          for(int k = 0; k < index.size(); k++){
            for(int l = 0; l < index.size(); l++){
              if(is_left[k] && is_right[l]){
                pairwise_tmrca[k][l] = coord;
                pairwise_tmrca[l][k] = coord;
              }
            }
          }
        }

      }


      //now determine Newick string of topologies of these three individuals
      int left_index, right_index, out_index;

      std::vector<std::string> newick(index.size(), "");
      for(int k = 0; k < index.size(); k++){
        newick[k] = sample.groups[sample.group_of_haplotype[index[k]]];
      }

      std::vector<bool> active(index.size(), true);
      int num_active = index.size();
      int c1 = 0, c2 = 1;
      std::vector<float> min(num_active-1, std::numeric_limits<float>::infinity());
      //float min = std::numeric_limits<float>::infinity();
      while(num_active >= 2){

        for(int i = 0; i < index.size(); i++){
          if(active[i]){
            for(int j = 0; j < index.size(); j++){
              if(active[j]){
                //std::cerr << pairwise_tmrca[i][j] << " " << min << std::endl; 
                if(i != j && pairwise_tmrca[i][j] < min[index.size() - num_active]){
                  c1 = i;
                  c2 = j;
                  min[index.size() - num_active] = pairwise_tmrca[i][j];
                }
              }
            }
          }
        }

        if(num_active == 3){
          left_index = index[c1];
          right_index = index[c2];
        }else if(num_active == 2){
          out_index = index[c2];
          if(out_index == left_index || out_index == right_index){
            out_index = index[c1];
          }
        }

        //std::cerr << c1 << " " << c2 << " " << min << std::endl;
        newick[c1] = "(" + newick[c1] + "," + newick[c2] + ")";
        active[c2] = false;
        num_active--;

      }	


      //record mutations
      //now need to check which branches are evidenced by a mutation
      //I have four types of branches (left, right, both, out)
      mut_left.clear();
      mut_right.clear();
      mut_both.clear();
      mut_out.clear();

      int ind1 = left_index, ind2 = right_index;
      while(1){
        while(ind2 < ind1){
          if(mut_on_branches[ind2].size() > 0){
            for(int k = 0; k < mut_on_branches[ind2].size(); k++){
              mut_right.push_back(mut_on_branches[ind2][k]);
            }
          }
          ind2 = (*mtr.tree.nodes[ind2].parent).label;
        }
        if(ind1 == ind2){
          break;
        }else{
          if(mut_on_branches[ind1].size() > 0){
            for(int k = 0; k < mut_on_branches[ind1].size(); k++){
              mut_left.push_back(mut_on_branches[ind1][k]);
            }
          }
          ind1 = (*mtr.tree.nodes[ind1].parent).label;
        }
      }
      assert(ind1 == ind2);

      ind2 = out_index;
      while(1){
        while(ind2 < ind1){
          if(mut_on_branches[ind2].size() > 0){
            for(int k = 0; k < mut_on_branches[ind2].size(); k++){
              mut_out.push_back(mut_on_branches[ind2][k]);
            }
          }
          ind2 = (*mtr.tree.nodes[ind2].parent).label;
        }
        if(ind1 == ind2){
          break;
        }else{
          if(mut_on_branches[ind1].size() > 0){
            for(int k = 0; k < mut_on_branches[ind1].size(); k++){
              mut_both.push_back(mut_on_branches[ind1][k]);
            }
          }
          ind1 = (*mtr.tree.nodes[ind1].parent).label;
        }
      }
      assert(ind1 == ind2);


      std::sort(mut_left.begin(), mut_left.end());
      std::sort(mut_right.begin(), mut_right.end());
      std::sort(mut_both.begin(), mut_both.end());
      std::sort(mut_out.begin(), mut_out.end());

      if(mut_left.size() + mut_right.size() + mut_both.size() + mut_out.size() > 0){
        os_loc_tmrca << chr << " " << (*it_mut).pos << " " << num_bases_tree_persists << " " << num_snps_on_tree << " " << frac_branches_with_snp << " " << newick[c1];
        for(int i = 0; i < min.size(); i++){
          os_loc_tmrca << " " << min[i];
        }
        os_loc_tmrca << "\n";

        for(int k = 0; k < mut_left.size(); k++){
          os_mut << chr << " " << (*it_mut).pos << " left " << ancmut.mut.info[mut_left[k]].pos << "\n";
        }
        for(int k = 0; k < mut_right.size(); k++){
          os_mut << chr << " " << (*it_mut).pos << " right " << ancmut.mut.info[mut_right[k]].pos << "\n";
        }
        for(int k = 0; k < mut_both.size(); k++){
          os_mut << chr << " " << (*it_mut).pos << " both " << ancmut.mut.info[mut_both[k]].pos << "\n";
        }
        for(int k = 0; k < mut_out.size(); k++){
          os_mut << chr << " " << (*it_mut).pos << " out " << ancmut.mut.info[mut_out[k]].pos << "\n";
        }
      }


      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      first_tree = false;
      count++;
    }

  }

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

void
TMRCAbyTopo2(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("haps") || !options.count("sample") || !options.count("input") || !options.count("output") || !options.count("poplabels")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output, poplabels." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate TMRCA by topology of pop_of_interest." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "TMRCA between groups ..." << std::endl;

  std::vector<std::string> filenames_anc, filenames_mut, filenames_chr, filenames_mask, filenames_haps, filenames_sample;

  if(options.count("chr")){
    igzstream is_chr(options["chr"].as<std::string>());
    if(is_chr.fail()){
      std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
    }
    std::string line;
    while(getline(is_chr, line)){
      filenames_anc.push_back(options["input"].as<std::string>() + "_chr" + line + ".anc"); 
      filenames_mut.push_back(options["input"].as<std::string>() + "_chr" + line + ".mut");    
      filenames_haps.push_back(options["haps"].as<std::string>() + "_chr" + line + ".haps.gz"); 
      filenames_sample.push_back(options["sample"].as<std::string>() + "_chr" + line + ".sample.gz"); 
      filenames_chr.push_back(line);
      if(options.count("mask")) filenames_mask.push_back(options["mask"].as<std::string>() + "_chr" + line + ".fa.gz");
    }
    is_chr.close();
  }else{
    filenames_anc.push_back(options["input"].as<std::string>() + ".anc"); 
    filenames_mut.push_back(options["input"].as<std::string>() + ".mut"); 
    filenames_haps.push_back(options["haps"].as<std::string>()); 
    filenames_sample.push_back(options["sample"].as<std::string>()); 
    filenames_chr.push_back("NA");  
    if(options.count("mask")) filenames_mask.push_back(options["mask"].as<std::string>());
  }

  Sample sample;
  sample.Read(options["poplabels"].as<std::string>());
  std::string label;

  if(!options.count("pop_of_interest")){
    label = sample.AssignPopOfInterest("All");
  }else{
    label = sample.AssignPopOfInterest(options["pop_of_interest"].as<std::string>());
  }
  assert(sample.group_of_interest.size() == 3);

  AncMutIterators ancmut_tmp(filenames_anc[0], filenames_mut[0]); 
  Data data(ancmut_tmp.NumTips(),ancmut_tmp.NumSnps());
  assert(data.N % 2 == 0);

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;

  std::vector<std::vector<int>> index(sample.group_of_interest.size());
  //choosing the first haplotype associated with group of interest
  for(int k = 0; k < sample.group_of_interest.size(); k++){
    std::cerr << sample.groups[sample.group_of_interest[k]] << std::endl;
    for(int i = 0; i < data.N; i++){
      if(sample.group_of_haplotype[i] == sample.group_of_interest[k]){
        index[k].push_back(i);
      }
    }
    std::sort(index[k].begin(), index[k].end());
  }

  ////////// 1. Read one tree at a time /////////

  //We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
  //The mut file is read once, file is closed after constructor is called.

  double genome_length = 0.0;
  std::vector<float> coordinates;
  std::vector<Leaves> desc;

  std::ofstream os_loc_tmrca(options["output"].as<std::string>() + ".topotmrca");
  if(os_loc_tmrca.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".topotmrca" << "." << std::endl;
    exit(1);
  }

  std::ofstream os_mut(options["output"].as<std::string>() + ".topomuts");
  if(os_mut.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".topomuts" << "." << std::endl;
    exit(1);
  }

  os_loc_tmrca << "CHR BP BP_persist num_snps_on_tree frac_branches_with_snp ind1 ind2 tmrca";
  for(int i = 0; i < index[0].size(); i++){
    os_loc_tmrca << " tmrca_" << sample.groups[sample.group_of_interest[1]] << "_with_" << i;
  }
  for(int i = 0; i < index[0].size(); i++){
    os_loc_tmrca << " tmrca_" << sample.groups[sample.group_of_interest[2]] << "_with_" << i;
  }
  os_loc_tmrca << "\n";

  std::vector<int> mut_left, mut_right, mut_both, mut_out;
  int count = 0;
  for(int chr = 0; chr < filenames_anc.size(); chr++){

    std::cerr << "chr: " << chr << std::endl;

    bool first_tree = true;

    //open genealogies
    AncMutIterators ancmut(filenames_anc[chr], filenames_mut[chr]); 
    assert(data.N == ancmut.NumTips());

    //open haps sample files
    haps m_hap(filenames_haps[chr].c_str(), filenames_sample[chr].c_str());
    if(data.N != ancmut.NumTips()){
      std::cerr << "Haps file and anc/mut have different number of samples" << std::endl;
      exit(1);
    }
    int num_snps = ancmut.NumSnps();

    /////////////////////////////////
    //Read anc and mut

    //read in the first 'len' SNPs in haps/sample
    int len = 1000000;
    std::vector<std::vector<char>> sequence(len);
    for(std::vector<std::vector<char>>::iterator it_seq = sequence.begin(); it_seq != sequence.end(); it_seq++){
      (*it_seq).resize(data.N);
    }
    std::vector<int> bp(len,0);
    int snp_mod = 0, snp_index = 0;
    for(snp_mod = 0; snp_mod < std::min(len, m_hap.GetL()); snp_mod++){
      m_hap.ReadSNP(sequence[snp_mod], bp[snp_mod]);
      snp_mod++;
    }
    snp_mod--;

    Leaves sequences_carrying_mutation;
    sequences_carrying_mutation.member.resize(data.N);

    /////////////////////////////////
    //Read anc and mut

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    while(num_bases_tree_persists >= 0.0){

      //new tree
      mtr.tree.GetCoordinates(coordinates);
      mtr.tree.FindAllLeaves(desc);

      //calculate some tree properties that correlate with tree quality
      float frac_branches_with_snp = 0, num_snps_on_tree = 0;
      for(std::vector<Node>::iterator it_n = std::next(mtr.tree.nodes.begin(), data.N); it_n != mtr.tree.nodes.end(); it_n++){
        num_snps_on_tree += (*it_n).num_events;
        if( ((*it_n).num_events) >= 1.0) frac_branches_with_snp += 1.0;
      }
      frac_branches_with_snp /= data.N-1.0;

      //by default, I know how many mutations map to each branch in the tree, but I'm not recording which specific mutations these are
      //So I need to remap these.

      //find the maximum range of the tree, i.e. each branch has a range in the genome, find maximum span
      int min_snp = mtr.tree.nodes[0].SNP_begin, max_snp = mtr.tree.nodes[0].SNP_end;
      for(int i = data.N; i < mtr.tree.nodes.size(); i++){
        if(min_snp > mtr.tree.nodes[i].SNP_begin) min_snp = mtr.tree.nodes[i].SNP_begin;
        if(max_snp < mtr.tree.nodes[i].SNP_end) max_snp = mtr.tree.nodes[i].SNP_end;
      }
      int min_bp = ancmut.mut.info[min_snp].pos, max_bp = ancmut.mut.info[max_snp].pos;

      //read all mutations between min_bp and max_bp
      int snp_start = -1, snp_end = -1, bp_dummy = bp[0];
      //bp[snp_mod] is last read entry, i.e. largest BP

      if(max_bp > bp[snp_mod]){ //need to parse additional sites in haps/sample
        while(max_bp > bp[snp_mod]){
          snp_mod++;
          snp_mod = snp_mod % len;
          m_hap.ReadSNP(sequence[snp_mod], bp[snp_mod]);
        }
      }

      //check if first entry is already exceeding
      if(bp[ (snp_mod + 1) % len ] >= min_bp){
        snp_start = (snp_mod + 1) % len;
        if(!first_tree) std::cerr << "Warnings: Haps/sample buffer may be too short." << std::endl;
      }else{
        snp_start = (snp_mod + 1) % len;
        for(int i = 0; i < len; i++){
          if(bp[ (i + snp_mod + 1) % len ] >= min_bp){
            snp_start = (i + snp_mod + 1) % len;
            break;
          }
        }
      }

      if(bp[ (snp_mod + 1) % len ] >= max_bp){
        snp_end = (snp_mod + 1) % len;
        std::cerr << "Warnings: Haps/sample buffer may be too short." << std::endl;
      }else{
        snp_end = (snp_mod + 1) % len;
        for(int i = 0; i < len; i++){
          if(bp[ (i + snp_mod + 1) % len ] >= max_bp){
            snp_end = (i + snp_mod + 1) % len;
            break;
          }
        }
      }

      //need to check for mutations between snp_start and snp_end
      //for each check whether they are mapping
      Tree tr = mtr.tree;
      std::vector<std::vector<int>> mut_on_branches(mtr.tree.nodes.size());
      int num_carriers;	

      //remap mutations
      //get snp index from (*it_mut)
      snp_index = 0;
      //map_mutations
      for(int snp = snp_start; snp <= snp_end; snp++){

        while(ancmut.mut.info[snp_index].pos < bp[snp]) snp_index++;
        if(ancmut.mut.info[snp_index].pos == bp[snp]){

          sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
          for(int i = 0; i < data.N; i++){
            if(sequence[snp % len][i] == '1'){
              sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
              sequences_carrying_mutation.num_leaves++;
            }else{
              sequences_carrying_mutation.member[i] = 0;
            }
          }

          if(sequences_carrying_mutation.num_leaves > 1){
            if(sequences_carrying_mutation.num_leaves > 0 && sequences_carrying_mutation.num_leaves < data.N){

              std::vector<SNPInfo> muts(1);
              Muts::iterator it_m = muts.begin();
              if(MapMutation(tr, sequences_carrying_mutation, it_m, false) < 3){
                int branch = (*it_m).branch[0];
                if(ancmut.mut.info[mtr.tree.nodes[branch].SNP_begin].pos <= bp[snp % len] && ancmut.mut.info[mtr.tree.nodes[branch].SNP_end].pos >= bp[snp % len]){
                  mut_on_branches[branch].push_back(snp_index);
                }
              }

            }
          }

          snp_index++;
        }

      }

      //remapping done

      std::vector<int> snplist1, snplist2;
      std::vector<float> tmrca1(index[0].size(), 0.0), tmrca2(index[0].size(), 0.0);

      for(int k1 = 0; k1 < index[1].size(); k1++){
        for(int k2 = 0; k2 < index[2].size(); k2++){
       
          //get TMRCA of index[1][k1] and index[2][k2]
          //also get a list of SNPs ancestral to each

          snplist1.clear();
          snplist2.clear();

          bool has_coal;
          float tmrca_12;
          int n, parent;
          int child;

          tmrca_12 = 0.0;
          has_coal = false;
          n = index[1][k1];
          while(tr.nodes[n].parent != NULL){
          
            if(!has_coal){
              for(int m = 0; m < mut_on_branches[n].size(); m++){
                snplist1.push_back(mut_on_branches[n][m]);
              }
            }

            parent = (*tr.nodes[n].parent).label;
            child = (*tr.nodes[parent].child_left).label;
            if(child == n){
              child = (*tr.nodes[parent].child_right).label;
            }

            for(std::vector<int>::iterator it_hap1 = desc[child].member.begin(); it_hap1 != desc[child].member.end(); it_hap1++){
              if(*it_hap1 == index[2][k2]){
                tmrca_12 = coordinates[parent];
                has_coal = true;
              }
              for(int k0 = 0; k0 < index[0].size(); k0++){
                if(*it_hap1 == index[0][k0]){
                  tmrca1[k0] = coordinates[parent]; 
                }
              }
            }
            n = parent;

          }

          tmrca_12 = 0.0;
          has_coal = false;
          n = index[2][k2];
          while(tr.nodes[n].parent != NULL){

            if(!has_coal){
              for(int m = 0; m < mut_on_branches[n].size(); m++){
                snplist2.push_back(mut_on_branches[n][m]);
              }
            }

            parent = (*tr.nodes[n].parent).label;
            child = (*tr.nodes[parent].child_left).label;
            if(child == n){
              child = (*tr.nodes[parent].child_right).label;
            }

            for(std::vector<int>::iterator it_hap1 = desc[child].member.begin(); it_hap1 != desc[child].member.end(); it_hap1++){
              if(*it_hap1 == index[1][k1]){
                tmrca_12 = coordinates[parent];
                has_coal = true;
              }
              for(int k0 = 0; k0 < index[0].size(); k0++){
                if(*it_hap1 == index[0][k0]){
                  tmrca2[k0] = coordinates[parent]; 
                }
              }
            }

            n = parent;

          }
          
          std::sort(tmrca1.begin(), tmrca1.end());
          std::sort(tmrca2.begin(), tmrca2.end());

          if(snplist1.size() > 0 && snplist2.size() > 0){

            //output  tmrca_12 tmrca1 tmrca2
            os_loc_tmrca << chr << " " << (*it_mut).pos << " " << num_bases_tree_persists << " " << num_snps_on_tree << " " << frac_branches_with_snp << " " << index[1][k1] << " " << index[2][k2] << " " << tmrca_12 << " ";
            for(std::vector<float>::iterator it = tmrca1.begin(); it != tmrca1.end(); it++){
              os_loc_tmrca << *it << " ";
            }
            for(std::vector<float>::iterator it = tmrca2.begin(); it != tmrca2.end(); it++){
              os_loc_tmrca << *it << " ";
            }
            os_loc_tmrca << "\n";

            for(std::vector<int>::iterator it1 = snplist1.begin(); it1 != snplist1.end(); it1++){
              os_mut << chr << " " << (*it_mut).pos << " " << sample.groups[sample.group_of_interest[1]] << " " << index[1][k1] << " " << ancmut.mut.info[*it1].pos << "\n";
            }
            for(std::vector<int>::iterator it2 = snplist2.begin(); it2 != snplist2.end(); it2++){
              os_mut << chr << " " << (*it_mut).pos << " " << sample.groups[sample.group_of_interest[2]] << " " << index[2][k2] << " " << ancmut.mut.info[*it2].pos << "\n";
            }

          }

        }
      }

      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      first_tree = false;
      count++;
    }

  }

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

