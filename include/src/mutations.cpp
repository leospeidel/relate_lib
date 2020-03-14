#include "mutations.hpp"


///////////////////////////////////////////////////

Mutations::Mutations(Data& data){
  N       = data.N;
  L       = data.L;

  info.resize(L);
  num_flips = 0; 
  num_notmappingmutations = 0; 
}

void
Mutations::Init(Data& data){
  N       = data.N;
  L       = data.L;

  info.resize(L);
  num_flips = 0; 
  num_notmappingmutations = 0; 
}

//////////////////// Public Members //////////

void
Mutations::GetAge(AncesTree& anc){

  CorrTrees::iterator it_seq = anc.seq.begin();
  int count_tree = 0;
  for(std::vector<SNPInfo>::iterator it = info.begin(); it != info.end(); it++){

    if((*it).tree > count_tree){
      count_tree++;
      it_seq++;
      assert((*it).tree == count_tree);
    }

    //find the branch on the tree and get age
    if((*it).branch.size() == 1){
      Node n = (*it_seq).tree.nodes[*(*it).branch.begin()];
      //traverse this node down to the bottom
      (*it).age_begin = 0.0;
      (*it).age_end = n.branch_length;
      while(n.child_left != NULL){
        n = *n.child_left;
        (*it).age_begin += n.branch_length;
      }
      (*it).age_end += (*it).age_begin;
    }
  }

}

void
Mutations::Read(igzstream& is){

  info.clear();
  info.resize(L);
  std::string line;
  int snp = 0;

  num_flips           = 0;
  num_notmappingmutations = 0;

  std::string inread;

  //file structure:
  //snp;pos_of_snp;rs-id;tree_index;branch_indices;is_mapping;is_flipped;(age_begin,age_end);ancestral_allele/alternative_allele;downstream_allele;upstream_allele;ACB;ASW;BEB;CDX;CEU;CHB;CHS;CLM;ESN;FIN;GBR;GIH;GWD;IBS;ITU;JPT;KHV;LWK;MSL;MXL;PEL;PJL;PUR;STU;TSI;YRI

  int i = 0;
  while(getline(is, line)){ 

    i = 0;

    //snp-id
    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    i++;
    info[snp].snp_id = std::stoi(inread);

    //pos_of_snp
    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    i++;
    info[snp].pos = std::stoi(inread);

    //dist_to_next_snp

    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    i++;
    info[snp].dist = std::stoi(inread);

    //rs-id
    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    i++;
    info[snp].rs_id = inread;

    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    info[snp].tree = std::stoi(inread);
    inread.clear();

    do{
      i++;
      while(line[i] != ' ' && line[i] != ';'){
        inread += line[i];
        i++;
      }
      if(inread.size() > 0){ 
        info[snp].branch.push_back(std::stoi(inread));
        inread.clear();
      }
    }while(line[i] != ';');

    i++;
    while(line[i] != ';') i++;
    i++;
    while(line[i] != ';'){
      inread += line[i];
      i++;
    } 
    info[snp].flipped = std::stoi(inread);
    inread.clear();
    i++;

    while(line[i] != ';'){
      inread += line[i];
      i++;
    } 
    info[snp].age_begin = std::stof(inread);
    inread.clear();
    i++;

    while(line[i] != ';'){
      inread += line[i];
      i++;
    } 
    info[snp].age_end = std::stof(inread);
    inread.clear();
    i++;

    if(i < line.size()){

      //allele type
      while(line[i] != ';' && i < line.size()){
        inread += line[i];
        i++;
      } 
      info[snp].mutation_type = inread;
      i++;
      if(i < line.size()){
        inread.clear();
        while(line[i] != ';'){
          inread += line[i];
          i++;
        } 
        info[snp].upstream_base = inread;
        i++;
        inread.clear();
        while(line[i] != ';'){
          inread += line[i];
          i++;
        } 
        info[snp].downstream_base = inread;
        i++;
        inread.clear();

        while(i < line.size()){ 
          inread.clear();
          while(line[i] != ';' && i < line.size()){
            inread += line[i];
            i++;
          } 
          i++;
          info[snp].freq.push_back(std::stoi(inread));
        }
      }

    }

    if(info[snp].flipped) num_flips++;
    if(info[snp].branch.size() > 1) num_notmappingmutations++;
    snp++;

  }

  info.resize(snp);
  is.close();

}

void
Mutations::Read(const std::string& filename){

  std::string line;
  igzstream is(filename);
  if(is.fail()) is.open(filename + ".gz");
  if(is.fail()){
    std::cerr << "Error while reading " << filename << "(.gz)." << std::endl;
    exit(1);
  }
  std::getline(is, header); 

  std::string unused;
  L = 0;
  while( std::getline(is, unused) ){
    ++L;
  }
  is.close();

  is.open(filename);
  if(is.fail()) is.open(filename + ".gz");
  std::getline(is, line);
  Read(is);

}

void 
Mutations::Dump(const std::string& filename){

  std::ofstream os;
  os.open(filename);

  if(os.fail()){
    std::cerr << "Error while writing to " << filename << "." << std::endl; 
  }else{

    if(header.size() > 0){
      os << header;
    }else{
      os << "snp;pos_of_snp;dist;rs-id;tree_index;branch_indices;is_not_mapping;is_flipped;age_begin;age_end;ancestral_allele/alternative_allele;upstream_allele;downstream_allele;"; 
    }
    os << "\n";

    for(std::vector<SNPInfo>::iterator it = info.begin(); it != info.end(); it++){
      os << (*it).snp_id << ";" << (*it).pos << ";" << (*it).dist << ";" << (*it).rs_id << ";" << (*it).tree << ";";
      std::deque<int>::iterator it_branch = (*it).branch.begin();
      if((*it).branch.size() > 0){
        os << *it_branch;
        it_branch++;
      }
      for(; it_branch != (*it).branch.end(); it_branch++){
        os << " " << *it_branch;
      }
      if((*it).branch.size() > 1){
        os << ";1;"; 
      }else{
        os << ";0;";
      }
      os << (*it).flipped << ";" << (*it).age_begin << ";" << (*it).age_end << ";";
      os << (*it).mutation_type << ";";

      if((*it).freq.size() > 0){
        os << (*it).upstream_base << ";" << (*it).downstream_base << ";";
        for(std::vector<int>::iterator it_freq = (*it).freq.begin(); it_freq != (*it).freq.end(); it_freq++){
          os << *it_freq << ";";
        }

      }

      os << "\n";

    }

    os.close();

  }

}

//////////////////////////////////////////

AncMutIterators::AncMutIterators(const std::string& filename_anc, const std::string& filename_mut): filename_anc(filename_anc), filename_mut(filename_mut){
  is.open(filename_anc);
  if(is.fail()) is.open(filename_anc + ".gz");
  if(is.fail()){
    std::cerr << "Failed to open file " << filename_anc << "(.gz)" << std::endl;
    exit(1);
  }

  //std::istringstream is_header, is_header2;

  std::string line, tmp;
  //read num_haplotypes
  getline(is, line);
  std::istringstream is_header(line);
  is_header >> tmp;
  is_header >> N;

  //read sample ages
  sample_ages.resize(N, 0.0);
  std::vector<double>::iterator it_sample_ages = sample_ages.begin();
  int i = 0;
  while(is_header >> *it_sample_ages){
    it_sample_ages++;
    i++;
    if(it_sample_ages == sample_ages.end()) break;
  }
  if(i != N) sample_ages.clear();

  //read num trees
  getline(is, line);
  std::istringstream is_header2(line);
  is_header2 >> tmp;
  is_header2 >> num_trees;

  mut.Read(filename_mut);
  pit_mut           = mut.info.begin();
  tree_index_in_mut = (*pit_mut).tree;
  tree_index_in_anc = -1;
}

void
AncMutIterators::OpenFiles(const std::string& i_filename_anc, const std::string& i_filename_mut){

  filename_anc = i_filename_anc;
  filename_mut = i_filename_mut;

  if(is.rdbuf() -> is_open()) is.close(); //close if stream is still open

  is.open(filename_anc);
  if(is.fail()) is.open(filename_anc + ".gz");
  if(is.fail()){
    std::cerr << "Failed to open file " << filename_anc << "(.gz)" << std::endl;
    exit(1);
  }

  std::istringstream is_header;

  std::string line, tmp;
  //read num_haplotypes
  getline(is, line);
  is_header.str(line);
  is_header >> tmp;
  is_header >> N;

  //read sample ages
  sample_ages.resize(N);
  std::vector<double>::iterator it_sample_ages = sample_ages.begin();
  int i = 0;
  while(is_header >> *it_sample_ages){
    it_sample_ages++;
    i++;
    if(it_sample_ages == sample_ages.end()) break;
  }
  if(i != N) sample_ages.clear();

  //read num trees
  getline(is, line);
  is_header.str(line);
  is_header.clear();
  is_header >> tmp;
  is_header >> num_trees;

  mut.Read(filename_mut);
  pit_mut           = mut.info.begin();
  tree_index_in_mut = (*pit_mut).tree;
  tree_index_in_anc = -1;

}

double
AncMutIterators::NextTree(MarginalTree& mtr, Muts::iterator& it_mut){

  if(tree_index_in_anc + 1 == num_trees){
    return(-1.0); //signals that we reached last tree
  }
  assert(getline(is, line));
  //read marginal tree
  if(sample_ages.size() > 0){
    mtr.Read(line, N, sample_ages);
  }else{
    mtr.Read(line, N);
  }
  tree_index_in_anc++;

  //pit_mut is pointing to first SNP in new tree already
  it_mut                         = pit_mut;

  if(tree_index_in_anc == tree_index_in_mut){

    //calculate how long tree persists 
    if(pit_mut != mut.info.begin()){
      assert((*std::prev(pit_mut,1)).dist >= 0.0);
      num_bases_tree_persists = (*std::prev(pit_mut,1)).dist/2.0;
    }else{
      num_bases_tree_persists = 0.0;
    }

    while(tree_index_in_mut == (*pit_mut).tree){
      assert((*pit_mut).dist >= 0.0);
      num_bases_tree_persists  += (*pit_mut).dist;
      pit_mut++;
      if(pit_mut == mut.info.end()) break;
    }

    if(pit_mut != mut.info.end()){
      assert((*std::prev(pit_mut,1)).dist >= 0.0);
      num_bases_tree_persists -= (*std::prev(pit_mut,1)).dist/2.0;   
      assert(tree_index_in_mut < (*pit_mut).tree);
      tree_index_in_mut        = (*pit_mut).tree;
    }

    return(num_bases_tree_persists);

  }else{ //tree has no mutations
    return 0.0;
  }

}

double
AncMutIterators::FirstSNP(MarginalTree& mtr, Muts::iterator& it_mut){

  if(tree_index_in_anc > 0){ //need to reopen file
    is.close();
    OpenFiles(filename_anc, filename_mut);
    NextTree(mtr, it_mut);
  }

  it_mut = mut.info.begin();

  if(it_mut == mut.info.end()){ //now at end
    is.close();
    return(-1.0);
  }

  //do nothing if (*it_mut).tree + 1 == tree_index_in_mut;
  if((*it_mut).tree == tree_index_in_mut){
    while(NextTree(mtr, it_mut) == 0.0); //skip any trees without mutations
  }

  //return the number of bases between prev and next snp (midpoints)
  if(it_mut != mut.info.begin()){
    assert((*std::prev(pit_mut,1)).dist >= 0.0);
    assert((*pit_mut).dist >= 0.0);
    num_bases_tree_persists = (*std::prev(it_mut,1)).dist/2.0 + (*it_mut).dist/2.0;
  }else{
    assert((*pit_mut).dist >= 0.0);
    num_bases_tree_persists = (*it_mut).dist/2.0;
  }

  return(num_bases_tree_persists);

}

double
AncMutIterators::NextSNP(MarginalTree& mtr, Muts::iterator& it_mut){

  if(it_mut == mut.info.end()){ //already at end
    is.close();
    return(-1.0);
  }
  it_mut++;
  if(it_mut == mut.info.end()){ //now at end
    is.close();
    return(-1.0);
  }

  //do nothing if (*it_mut).tree + 1 == tree_index_in_mut;
  if((*it_mut).tree == tree_index_in_mut){
    while(NextTree(mtr, it_mut) == 0.0); //skip any trees without mutations
  }

  //return the number of bases between prev and next snp (midpoints)
  if(it_mut != mut.info.begin()){
    assert((*std::prev(pit_mut,1)).dist >= 0.0);
    assert((*pit_mut).dist >= 0.0);
    num_bases_tree_persists = (*std::prev(it_mut,1)).dist/2.0 + (*it_mut).dist/2.0;
  }else{
    assert((*pit_mut).dist >= 0.0);
    num_bases_tree_persists = (*it_mut).dist/2.0;
  }

  return(num_bases_tree_persists);

}


