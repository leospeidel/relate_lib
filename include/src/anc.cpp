#include "anc.hpp"

//////////////////////////////////

void
Tree::ReadTree(const char* line, int N){

  //clock_t tbegin = clock();
  int i = 0;
  nodes.clear();
  nodes.resize(2*N-1);

  std::vector<Node>::iterator it_node = nodes.begin(); 
  Node* p_parent;
  int num_node = 0, parent;

  while(num_node < 2*N-1){

    sscanf(&line[i], "%d:(%lf %f %d %d)", &parent, &(*it_node).branch_length, &(*it_node).num_events, &(*it_node).SNP_begin, &(*it_node).SNP_end); 

    while(line[i] != ')' && line[i] != '\n') i++;
    i += 2;

    if(parent != -1){
      p_parent   = &nodes[parent];
      (*it_node).parent         = p_parent;
      (*it_node).label          = num_node;
      if((*p_parent).child_left == NULL){
        (*p_parent).child_left  = &(*it_node);
      }else{
        (*p_parent).child_right = &(*it_node);       
      }
    }else{
      (*it_node).parent         = NULL;
      (*it_node).label          = num_node;
    }       

    num_node++;
    it_node++;

  }

  //clock_t tend = clock();
  //double elapsed_secs = double(tend - tbegin) / CLOCKS_PER_SEC;
  //std::cerr << elapsed_secs << std::endl;

}

void
Tree::WriteNewick(const std::string& filename_newick, double factor, const bool add) const{

  //coordinates.clear();

  //maybe not the most efficient convertion algorithm but it works

  int root = (int) nodes.size() - 1;
  for(int i = 0; i < (int) nodes.size(); i++){
    if(nodes[i].parent == NULL){
      root = i; //root
      break;
    }
  }

  std::ofstream os_new;
  if(!add){ 
    os_new.open(filename_newick);
  }else{
    os_new.open(filename_newick, std::ofstream::app);
  }

  std::list<Node> todo_nodes;
  float l1 = ((*nodes[root].child_left).branch_length) * factor;
  float l2 = ((*nodes[root].child_right).branch_length) * factor;

  //if(l1 < 0.0) std::cerr << root << ": " << coordinates[root].tree << ", " << children[root][0] << ": " << coordinates[children[root][0]].tree << std::endl;    
  //if(l2 < 0.0) std::cerr << root << ": " << coordinates[root].tree << ", " << children[root][1] << ": " << coordinates[children[root][1]].tree << std::endl;

  std::string newick = "(" + std::to_string((*nodes[root].child_left).label) + ":" + std::to_string(l1)  + "," + std::to_string((*nodes[root].child_right).label) + ":" + std::to_string(l2) + ")";
  todo_nodes.push_back(*nodes[root].child_left);
  todo_nodes.push_back(*nodes[root].child_right);

  while(todo_nodes.size() > 0){

    std::list<Node>::iterator node = todo_nodes.begin();

    if((*node).child_left == NULL){
      todo_nodes.erase(node);
    }else{
      int index = 0;
      std::string node_index;
      for(; index < (int) newick.size(); index++){
        while( !isdigit(newick[index]) ){
          //std::cerr << newick[index] << " " << isdigit(newick[index]) << std::endl;
          index++;
        }
        node_index.clear();
        while( newick[index] != ',' && newick[index] != ')' && newick[index] != ':'){
          node_index = node_index + newick[index];
          index++;
          //std::cerr << node_index << " " << index << " " << newick[index] << std::endl;
        }
        if((float)(*node).label == stof(node_index) && newick[index] == ':') break;
      }
      assert(index < (int) newick.size());

      float l1 = ((*(*node).child_left).branch_length) * factor;
      float l2 = ((*(*node).child_right).branch_length) * factor;

      std::string new_brackets = "(" + std::to_string((*(*node).child_left).label) + ":" + std::to_string(l1) + "," + std::to_string((*(*node).child_right).label) + ":" + std::to_string(l2) + ")";
      newick.replace(index-node_index.size(), node_index.size(), new_brackets);
      todo_nodes.push_back(*(*node).child_left);
      todo_nodes.push_back(*(*node).child_right);
      todo_nodes.erase(node);

    }

    //std::cerr << "newick: " << newick << " " << todo_nodes.size() << std::endl;

  }
  newick += ";";

  os_new << newick << std::endl;

  os_new.close();


}

void
Tree::WriteNHX(const std::string& filename_nhx, std::vector<std::string>& property, const bool add) const{

  if(property.size() != nodes.size()){
    std::cerr << "Property vector has wrong size." << std::endl;
    exit(1);
  }

  //maybe not the most efficient convertion algorithm but it works
  int root = (int) nodes.size() - 1;
  for(int i = 0; i < (int) nodes.size(); i++){
    if(nodes[i].parent == NULL){
      root = i; //root
      break;
    }
  }

  std::ofstream os_new;
  if(!add){ 
    os_new.open(filename_nhx);
  }else{
    os_new.open(filename_nhx, std::ofstream::app);
  }

  std::list<Node> todo_nodes; //stores node indices in newick string
  float l1 = ((*nodes[root].child_left).branch_length);
  float l2 = ((*nodes[root].child_right).branch_length);

  std::string newick = "(" + std::to_string((*nodes[root].child_left).label) + ":" + std::to_string(l1) + "[&&NHX:S=" + property[(*nodes[root].child_left).label] + "]," + std::to_string((*nodes[root].child_right).label) + ":" + std::to_string(l2) + "[&&NHX:S=" + property[(*nodes[root].child_right).label] + "])";
  todo_nodes.push_back(*nodes[root].child_left);
  todo_nodes.push_back(*nodes[root].child_right);

  while(todo_nodes.size() > 0){

    std::list<Node>::iterator node = todo_nodes.begin();

    if((*node).child_left == NULL){
      todo_nodes.erase(node); //erase if it is a tip
    }else{
      int index = 0;
      std::string node_index;
      for(; index < (int) newick.size(); index++){
        while( !isdigit(newick[index]) ){
          index++;
        }
        node_index.clear();
        while( newick[index] != ',' && newick[index] != ')' && newick[index] != ':'){
          node_index = node_index + newick[index];
          index++;
        }
        if((float)(*node).label == stof(node_index) && newick[index] == ':') break;
      }
      assert(index < (int) newick.size());

      float l1 = ((*(*node).child_left).branch_length);
      float l2 = ((*(*node).child_right).branch_length);

      std::string new_brackets = "(" + std::to_string((*(*node).child_left).label) + ":" + std::to_string(l1) + "[&&NHX:S=" + property[(*(*node).child_left).label]  + "]," + std::to_string((*(*node).child_right).label) + ":" + std::to_string(l2) + "[&&NHX:S=" + property[(*(*node).child_right).label] + "])";
      newick.replace(index-node_index.size(), node_index.size(), new_brackets); //replace node index by branching event
      todo_nodes.push_back(*(*node).child_left);
      todo_nodes.push_back(*(*node).child_right);
      todo_nodes.erase(node);

    }

  }
  newick += "[&&NHX:S=" + property[root] + "];";

  os_new << newick << std::endl;

  os_new.close();

}

////////////
void 
Tree::FindAllLeaves(std::vector<Leaves>& leaves) const{

  int N_total = nodes.size();
  int N = (N_total + 1)/2;
  leaves.resize(N_total);

  Node root = nodes[N_total - 1];
  if(root.parent != NULL){
    for(int i = N; i < N_total; i++){
      if(nodes[i].parent == NULL){
        root = nodes[i];
        break;
      }
    }
  }

  FindLeaves(root, leaves);

}

void 
Tree::FindLeaves(Node& node, std::vector<Leaves>& leaves) const{

  if(node.child_left != NULL){

    Node child1 = *node.child_left;
    Node child2 = *node.child_right;

    FindLeaves(child1, leaves);
    FindLeaves(child2, leaves);

    leaves[node.label].member.resize( leaves[child1.label].member.size() + leaves[child2.label].member.size() );

    std::vector<int>::iterator it_member = leaves[node.label].member.begin();
    std::vector<int>::iterator it_child1_member = leaves[child1.label].member.begin();
    std::vector<int>::iterator it_child2_member = leaves[child2.label].member.begin();

    const std::vector<int>::iterator it_child1_member_end = leaves[child1.label].member.end();
    const std::vector<int>::iterator it_child2_member_end = leaves[child2.label].member.end();

    for(; it_member != leaves[node.label].member.end();){

      if(it_child1_member != it_child1_member_end && it_child2_member != it_child2_member_end){
        if(*it_child1_member < *it_child2_member){
          *it_member = *it_child1_member;
          it_child1_member++;
          it_member++;
        }else{
          *it_member = *it_child2_member;
          it_child2_member++;
          it_member++;
        }
      }else if(it_child1_member != it_child1_member_end){
        *it_member = *it_child1_member;
        it_child1_member++;
        it_member++;
      }else{
        *it_member = *it_child2_member;
        it_child2_member++;
        it_member++;
      }

    }

    leaves[node.label].num_leaves = leaves[child1.label].num_leaves + leaves[child2.label].num_leaves;

  }else{
    leaves[node.label].member.resize(1);    //resizing bitset and filling it with FALSE
    leaves[node.label].member[0]  = node.label; //setting position node.label to TRUE
    leaves[node.label].num_leaves = 1;
  }

}

void
Tree::TraverseTreeToGetCoordinates(Node& n, std::vector<float>& coordinates){

  if(n.child_left != NULL){

    TraverseTreeToGetCoordinates(*n.child_left, coordinates);
    TraverseTreeToGetCoordinates(*n.child_right, coordinates);
    coordinates[n.label] = coordinates[(*n.child_left).label] + (*n.child_left).branch_length;  

    if(coordinates[n.label] < coordinates[(*n.child_right).label]){
      coordinates[n.label] = coordinates[(*n.child_right).label];
    }

  }else{
    coordinates[n.label] = 0.0;
  }

}

void 
Tree::GetCoordinates(std::vector<float>& coordinates){
  coordinates.resize(nodes.size());
  TraverseTreeToGetCoordinates(nodes[nodes.size() - 1], coordinates);
}

void 
Tree::GetCoordinates(int node, std::vector<float>& coordinates){
  TraverseTreeToGetCoordinates(nodes[node], coordinates);
}

/////////////////////////////////

void
Tree::GetNumberOfLeavesInSubpop(const Node& n, std::vector<int>& subpop, std::vector<int>& number_in_subpop) const{

  if(n.child_left != NULL){

    GetNumberOfLeavesInSubpop((*n.child_left), subpop, number_in_subpop);
    GetNumberOfLeavesInSubpop((*n.child_right), subpop, number_in_subpop);
    number_in_subpop[n.label] = number_in_subpop[(*n.child_left).label] + number_in_subpop[(*n.child_right).label];

  }else{

    for(std::vector<int>::iterator it_subpop = subpop.begin(); it_subpop != subpop.end(); it_subpop++){
      if(*it_subpop == n.label){
        number_in_subpop[*it_subpop] = 1;
        break;
      }
    }

  }

}

void
Tree::GetSubTree(Sample& sample, Tree& subtree) const{

  std::vector<int> subpop;

  int hap = 0;
  for(std::vector<int>::iterator it_group_of_haplotype = sample.group_of_haplotype.begin(); it_group_of_haplotype != sample.group_of_haplotype.end(); it_group_of_haplotype++){
    bool exists = false;
    for(std::vector<int>::iterator it_group_of_interest = sample.group_of_interest.begin(); it_group_of_interest != sample.group_of_interest.end(); it_group_of_interest++){
      if(*it_group_of_haplotype == *it_group_of_interest){
        exists = true;
        break;
      }
    }
    if(exists){
      subpop.push_back(hap);
    }
    hap++;
  }

  std::vector<int> convert_index, number_in_subpop;
  GetSubTree(subpop, subtree, convert_index, number_in_subpop);

}

void
Tree::GetSubTree(Sample& sample, Tree& subtree,  std::vector<int>& convert_index, std::vector<int>& number_in_subpop) const{

  std::vector<int> subpop;

  int hap = 0;
  for(std::vector<int>::iterator it_group_of_haplotype = sample.group_of_haplotype.begin(); it_group_of_haplotype != sample.group_of_haplotype.end(); it_group_of_haplotype++){
    bool exists = false;
    for(std::vector<int>::iterator it_group_of_interest = sample.group_of_interest.begin(); it_group_of_interest != sample.group_of_interest.end(); it_group_of_interest++){
      if(*it_group_of_haplotype == *it_group_of_interest){
        exists = true;
        break;
      }
    }
    if(exists){
      subpop.push_back(hap);
    }
    hap++;
  }

  GetSubTree(subpop, subtree, convert_index, number_in_subpop);

}

void
Tree::GetSubTree(std::vector<int>& subpop, Tree& subtree, std::vector<int>& convert_index, std::vector<int>& number_in_subpop) const{

  convert_index.resize(nodes.size(), -1);
  std::fill(convert_index.begin(), convert_index.end(), -1);

  //for each node, calculate the number of leaves in subpop below it
  number_in_subpop.resize(nodes.size(), 0.0);
  std::fill(number_in_subpop.begin(), number_in_subpop.end(), 0.0);
  GetNumberOfLeavesInSubpop(*std::prev(nodes.end(),1), subpop, number_in_subpop);

  if(subpop.size() >= (int)(nodes.size() + 1.0)/2.0){

    subtree = *this;
    for(int i = 0; i < (int)nodes.size(); i++){
      convert_index[i] = i;
    }

  }else{

    int node = 0;
    subtree.nodes.resize(2*subpop.size() - 1);

    for(node = 0; node < subpop.size(); node++){
      subtree.nodes[node] = nodes[subpop[node]];
      subtree.nodes[node].label = node;
      convert_index[subpop[node]] = node;
    }

    for(int i = (int)(nodes.size() + 1.0)/2.0; i < (int) nodes.size(); i++){

      int child_left  = (*nodes[i].child_left).label;
      int child_right = (*nodes[i].child_right).label;

      if(number_in_subpop[child_left] > 0 && number_in_subpop[child_right] > 0){

        assert(convert_index[child_left] != -1);
        assert(convert_index[child_right] != -1);
        subtree.nodes[node] = nodes[i];
        //connect to children.
        subtree.nodes[node].label       = node;
        subtree.nodes[node].child_left  = &subtree.nodes[convert_index[child_left]];
        subtree.nodes[node].child_right = &subtree.nodes[convert_index[child_right]];
        subtree.nodes[convert_index[child_left]].parent  = &subtree.nodes[node];
        subtree.nodes[convert_index[child_right]].parent = &subtree.nodes[node];

        convert_index[i]    = node;
        assert(number_in_subpop[i] > 0);
        node++;

      }else if(number_in_subpop[child_left] > 0){

        assert(convert_index[child_left] != -1);
        convert_index[i]                               = convert_index[child_left];
        subtree.nodes[convert_index[i]].branch_length += nodes[i].branch_length;
        subtree.nodes[convert_index[i]].num_events    += nodes[i].num_events;
        assert(number_in_subpop[i] > 0);

      }else if(number_in_subpop[child_right] > 0){

        assert(convert_index[child_right] != -1);
        convert_index[i]                               = convert_index[child_right];
        subtree.nodes[convert_index[i]].branch_length += nodes[i].branch_length;
        subtree.nodes[convert_index[i]].num_events    += nodes[i].num_events;
        assert(number_in_subpop[i] > 0);

      }else{

        assert(convert_index[i] == -1);
        assert(number_in_subpop[i] == 0);

      }

    }

    assert(node == 2*subpop.size()-1);
    subtree.nodes[node-1].parent = NULL;
  }

}


/////////////////////////////////

void
MarginalTree::Read(const std::string& line, int N){

  int i = 0;
  while(line[i] != ':'){
    i++;
  }
  i += 2;

  sscanf(line.c_str(), "%d: ", &pos);
  tree.ReadTree(&(line.c_str())[i], N);

}


////////////////////////////////

void 
AncesTree::Read(igzstream& is){

  std::string line;
  seq.clear();
  seq.emplace_back();
  CorrTrees::iterator it_seq = seq.begin();

  int num_tree = 0;
  int i;
  while(num_tree < L){

    getline(is, line);

    i = 0;
    while(line[i] != ':'){
      i++;
    }
    i += 2;
    sscanf(line.c_str(), "%d: ", &(*it_seq).pos);
    (*it_seq).tree.ReadTree(&(line.c_str())[i], N);

    seq.emplace_back();
    it_seq++;
    num_tree++;
  }

  seq.pop_back();

}

void
AncesTree::Read(const std::string& filename){

  double start_time = time(NULL);
  clock_t begin = clock();

  igzstream is(filename);
  if(is.fail()) is.open(filename + ".gz");
  if(is.fail()){ 
    std::cerr << "Error while opening file " << filename << "(.gz)." << std::endl;
    exit(1);
  }

  char sdummy[30];
  std::string line;
  getline(is, line);
  sscanf(line.c_str(), "%s %d", sdummy, &N);
  getline(is, line);
  sscanf(line.c_str(), "%s %d", sdummy, &L);

  Read(is);

  clock_t end = clock();
  double end_time = time(NULL);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  //std::cerr << "anc read in: " << elapsed_secs << " CPU secs and " << end_time - start_time << " real secs." << std::endl;

}


void
AncesTree::Dump(FILE *pfile){

  double start_time = time(NULL);
  clock_t begin = clock();

  int parent;
  std::vector<Node>::iterator n_it;
  for(CorrTrees::iterator it_seq = seq.begin(); it_seq != seq.end(); it_seq++){

    n_it = (*it_seq).tree.nodes.begin();
    fprintf(pfile, "%d: ", (*it_seq).pos);
    for(; n_it != (*it_seq).tree.nodes.end(); n_it++){
      if((*n_it).parent == NULL){
        parent = -1;
      }else{
        parent = (*(*n_it).parent).label;
      }
      fprintf(pfile, "%d:(%.5f %.3f %d %d) ", parent, (*n_it).branch_length, (*n_it).num_events, (*n_it).SNP_begin, (*n_it).SNP_end);
    }

    fprintf(pfile, "\n");

  }

  clock_t end = clock();
  double end_time = time(NULL);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  //std::cerr << "anc dumped in: " << elapsed_secs  << " CPU secs and " << end_time - start_time << " real secs." << std::endl;

}

void
AncesTree::Dump(std::ofstream& os){

  int parent;
  std::vector<Node>::iterator n_it;
  //std::string line;
  std::stringstream stream;
  for(CorrTrees::iterator it_seq = seq.begin(); it_seq != seq.end(); it_seq++){

    n_it = (*it_seq).tree.nodes.begin();
    stream << (*it_seq).pos << ": ";
    for(; n_it != (*it_seq).tree.nodes.end(); n_it++){
      if((*n_it).parent == NULL){
        parent = -1;
      }else{
        parent = (*(*n_it).parent).label;
      }
      stream << parent << ":(" << std::setprecision(5) << (*n_it).branch_length << " "  << std::setprecision(2) << (*n_it).num_events << " " << (*n_it).SNP_begin << " " << (*n_it).SNP_end << ")";
    }

    stream << "\n";
  }

  os << stream.str();

}

void
AncesTree::Dump(const std::string& filename){

  FILE *pfile = std::fopen(filename.c_str(), "w");

  if(pfile == NULL){

    std::cerr << "Error while writing to " << filename << "." << std::endl;

  }else{

    fprintf(pfile, "NUM_HAPLOTYPES %ld\n", ((*seq.begin()).tree.nodes.size() + 1)/2);
    fprintf(pfile, "NUM_TREES %ld\n", seq.size());

    Dump(pfile);

    fclose(pfile);

  }

}


////////////////////////////////////

void
AncesTree::ReadArgweaverSMC(const std::string& filename){

  int N       = 0;
  int N_total = 0;

  igzstream is(filename);
  std::string line;

  std::string sdummy, newick;
  int idummy;
  int node_index;
  std::vector<int> convert_index;

  //get number of nodes
  std::getline(is, line);
  std::stringstream readN(line);
  readN >> sdummy;
  while(readN >> idummy){
    convert_index.push_back(idummy-1);
    N++;
  }
 
  N_total = 2*N-1;
  convert_index.resize(N_total);
  for(int n = N; n < N_total;n++){
    convert_index[n] = n;
  }
  //std::cerr << N << std::endl;

  //std::getline(is, line);
  int num_tree = 0;
  seq.emplace_back();
  CorrTrees::iterator it_seq = seq.begin();
  while(std::getline(is, line)){ 

    assert(std::getline(is, line));

    std::stringstream tree_stream(line); 
    (*it_seq).tree.nodes.resize(N_total);

    tree_stream >> sdummy;
    tree_stream >> (*it_seq).pos;
    tree_stream >> idummy;

    //std::cerr << (*it_seq).pos << std::endl; 
    //(*it_seq).pos -= 1.0;
    tree_stream >> newick; 

    //need to convert newick into my tree data structure
    int i = 0;
    while(newick.size() > 0){ 
      int startpos, endpos;
      std::string index_child1, index_child2, index_parent;
      std::string branch_length1, branch_length2;

      //std::cerr << newick << std::endl << std::endl;
      while(newick[i] == '(') i++;
      startpos = i-1;
      while(newick[i] != ':'){
        index_child1 += newick[i];
        i++;
      }
      i++;
      while(newick[i] != '['){
        branch_length1 += newick[i];
        i++;
      }
      while(newick[i] != ',') i++;
      i++;
      if(newick[i] != '('){
        while(newick[i] != ':'){
          index_child2 += newick[i];
          i++;
        } 
        i++;
        while(newick[i] != '['){
          branch_length2 += newick[i];
          i++;
        }
        while(newick[i] != ')') i++;
        i++;
        endpos = i;
        while( newick[i] != ':' && newick[i] != '['){
          index_parent += newick[i];
          i++;
        }
        //std::cerr << newick[startpos] << " " << newick[endpos-1] << " " << newick[i] << std::endl; 
        //std::cerr << index_child1 << " " << branch_length1 << " " << index_child2 << " " << branch_length2 << " " << index_parent << std::endl;

        int child_left = convert_index[stoi(index_child1)], child_right = convert_index[stoi(index_child2)], parent = convert_index[stoi(index_parent)];
        //std::cerr << child_left << " " << bl1 << " " << child_right << " " << bl2 << " " << parent << std::endl << std::endl;

        (*it_seq).tree.nodes[child_left].label  = child_left;
        (*it_seq).tree.nodes[child_right].label = child_right;
        (*it_seq).tree.nodes[parent].label      = parent;

        (*it_seq).tree.nodes[child_left].parent  = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[child_right].parent = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[child_left];
        (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[child_right];

        (*it_seq).tree.nodes[child_left].branch_length = stof(branch_length1);
        (*it_seq).tree.nodes[child_right].branch_length = stof(branch_length2);

        if(newick[i] == '['){ 
          break;
        }
        newick.replace(startpos, endpos - startpos, "");  
        i = 0;
      }

    }

    //make sure that node root is root
    //find root
    int root = N_total-1;
    if((*it_seq).tree.nodes[root].parent != NULL){

      int real_root = 0;
      for(std::vector<Node>::iterator it_node = (*it_seq).tree.nodes.begin(); it_node != (*it_seq).tree.nodes.end(); it_node++){
        if((*it_node).parent == NULL){
          break;
        }
        real_root++;
      }

      Node n1 = (*it_seq).tree.nodes[real_root];
      Node n2 = (*it_seq).tree.nodes[root];

      bool left = false;
      if((*(*n2.parent).child_left).label == root){
        left = true;
      }else{ 
        assert((*(*n2.parent).child_right).label == root);
      }

      if((*n2.parent).label == real_root){

        (*it_seq).tree.nodes[real_root]       = n2;
        (*it_seq).tree.nodes[real_root].label = real_root; 
        (*it_seq).tree.nodes[root]       = n1;
        (*it_seq).tree.nodes[root].label = root;

        if(left){
          (*it_seq).tree.nodes[root].child_left  = &(*it_seq).tree.nodes[real_root];
        }else{ 
          (*it_seq).tree.nodes[root].child_right = &(*it_seq).tree.nodes[real_root];
        }
        (*(*it_seq).tree.nodes[real_root].child_left).parent  = &(*it_seq).tree.nodes[real_root];
        (*(*it_seq).tree.nodes[real_root].child_right).parent = &(*it_seq).tree.nodes[real_root];
        (*(*it_seq).tree.nodes[root].child_left).parent       = &(*it_seq).tree.nodes[root];
        (*(*it_seq).tree.nodes[root].child_right).parent      = &(*it_seq).tree.nodes[root];

      }else{

        (*it_seq).tree.nodes[real_root]       = n2;
        (*it_seq).tree.nodes[real_root].label = real_root; 
        if(left){
          (*(*it_seq).tree.nodes[real_root].parent).child_left  = &(*it_seq).tree.nodes[real_root];
        }else{ 
          (*(*it_seq).tree.nodes[real_root].parent).child_right = &(*it_seq).tree.nodes[real_root];
        }
        (*(*it_seq).tree.nodes[real_root].child_left).parent  = &(*it_seq).tree.nodes[real_root];
        (*(*it_seq).tree.nodes[real_root].child_right).parent = &(*it_seq).tree.nodes[real_root];

        (*it_seq).tree.nodes[root]       = n1;
        (*it_seq).tree.nodes[root].label = root;
        (*(*it_seq).tree.nodes[root].child_left).parent  = &(*it_seq).tree.nodes[root];
        (*(*it_seq).tree.nodes[root].child_right).parent = &(*it_seq).tree.nodes[root];
   
      }
      
      /* 
      std::cerr << (*it_seq).tree.nodes[root].label << " " << ((*it_seq).tree.nodes[root].parent == NULL) << " " << (*(*(*it_seq).tree.nodes[root].child_left).parent).label << " " << (*(*(*it_seq).tree.nodes[root].child_right).parent).label << std::endl;
  
      std::cerr << (*it_seq).tree.nodes[real_root].label << " " << ((*it_seq).tree.nodes[real_root].parent == NULL) << " " << (*(*(*it_seq).tree.nodes[real_root].child_left).parent).label << " " << (*(*(*it_seq).tree.nodes[real_root].child_right).parent).label << std::endl;

      std::cerr << (*(*(*it_seq).tree.nodes[real_root].parent).child_left).label << " " << (*(*(*it_seq).tree.nodes[real_root].parent).child_right).label << std::endl;

      std::cerr << "------------" << std::endl;
      if(num_tree == 789) exit(1);     
      */
      

    }




    num_tree++;
    seq.emplace_back();
    it_seq++;
    //std::getline(is, line); //skipping every second line
  }

  seq.pop_back();

  //std::cerr << seq.size() << std::endl;
}

void
AncesTree::ReadRent(const std::string& filename, float Ne){

  igzstream is(filename);
  std::string line;

  std::string newick; 

  int N = -1;
  int N_total;
  int i;

  std::vector<float> coordinates;
  int num_tree = 1;
  seq.resize(1);
  CorrTrees::iterator it_seq = seq.begin();
  while(std::getline(is, line)){ 

    i = 0;
    if(N == -1){
      N = 0;
      while(i < line.size()){
        if(line[i] == ',') N++;
        i++;
      }
      N += 1;
      N_total = 2*N-1;

    }

    //std::cerr << line << " " << N_total << std::endl;
    //exit(1);

    std::stringstream tree_stream(line); 
    (*it_seq).tree.nodes.resize(N_total);
    tree_stream >> (*it_seq).pos;
    tree_stream >> newick; 

    //need to convert newick into my tree data structure
    i = 0;
    int node = N+1;
    int count_bracket = 0, count_comma = 0;
    while(node <= N_total){ 
      int startpos, endpos;
      std::string index_child1, index_child2, index_parent;
      std::string branch_length1, branch_length2;

      while(newick[i] == '(') i++;
      startpos = i;

      while(newick[i] != ':'){
        index_child1 += newick[i];
        i++;
      }
      i++;
      while(newick[i] != ',' && i < newick.size()){
        branch_length1 += newick[i];
        i++;
      }
      i++;
      if(newick[i] != '(' && i < newick.size()){
        while(newick[i] != ':'){
          index_child2 += newick[i];
          i++;
        } 
        i++;
        while(newick[i] != ')' && i < newick.size()){
          branch_length2 += newick[i];
          i++;
        }
        i++;
        endpos = i;

        //std::cerr << index_child1 << " " << index_child2 << " " << branch_length1 << " " << branch_length2 << std::endl;
        int child_left = stoi(index_child1)-1, child_right = stoi(index_child2)-1, parent = node-1;

        (*it_seq).tree.nodes[child_left].label  = child_left;
        (*it_seq).tree.nodes[child_right].label = child_right;
        (*it_seq).tree.nodes[parent].label      = parent;

        (*it_seq).tree.nodes[child_left].parent  = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[child_right].parent = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[child_left];
        (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[child_right];

        (*it_seq).tree.nodes[child_left].branch_length  = stof(branch_length1) * Ne;
        (*it_seq).tree.nodes[child_right].branch_length = stof(branch_length2) * Ne;

        //std::cerr << newick << std::endl;
        //std::cerr << newick.substr(startpos-1, endpos-startpos+1) << std::endl;
        newick.replace(startpos-1, endpos - startpos + 1, std::to_string(node)); 
        count_bracket = 0; 
        count_comma = 0;
        i = 0;
        for(; i < newick.size(); i++){
          if(newick[i] == '(') count_bracket++;
          if(newick[i] == ',') count_comma++;
        }
        //std::cerr << count_comma << " " << count_bracket << std::endl;
        //assert(count_comma == count_bracket);
        if(count_comma != count_bracket){
          break;
        }
        //std::cerr << newick << std::endl; 
        i = 0;
        node++;
      }
    }

    bool everyone_has_parent = true;
    for(int i = 0; i < N_total - 1; i++){
      if((*it_seq).tree.nodes[i].parent == NULL){
        everyone_has_parent = false;
        break;
      } 
    }
    //std::cerr << everyone_has_parent << std::endl;
    if(node != N_total + 1 || count_comma != count_bracket || !everyone_has_parent){
      seq.pop_back();
      num_tree--;
      it_seq--;
    }

    /*
       coordinates.resize(N_total); 
       (*it_seq).tree.GetCoordinates(coordinates);
       std::cerr << coordinates[N_total-1] << " " << coordinates[N_total-1]/Ne << std::endl; 
       */

    num_tree++;
    seq.emplace_back();
    it_seq = std::prev(seq.end(),1);
    //it_seq++;
  }

  seq.pop_back();

  //std::cerr << num_tree << " " << (*std::prev(seq.end(),1)).pos << std::endl;
}

void
AncesTree::ReadNewick(const std::string& filename, float Ne){

  igzstream is(filename);
  std::string line;

  std::string newick; 

  int N = -1;
  int N_total;
  int i;

  std::vector<float> coordinates;
  int num_tree = 1;
  seq.resize(1);
  CorrTrees::iterator it_seq = seq.begin();
  while(std::getline(is, line)){ 

    i = 0;
    if(N == -1){
      N = 0;
      while(i < line.size()){
        if(line[i] == ',') N++;
        i++;
      }
      N += 1;
      N_total = 2*N-1;

    }

    //std::cerr << line << " " << N_total << std::endl;
    //exit(1);

    std::stringstream tree_stream(line); 
    (*it_seq).tree.nodes.resize(N_total);
    tree_stream >> (*it_seq).pos;
    tree_stream >> newick; 

    //need to convert newick into my tree data structure
    i = 0;
    int node = N;
    int count_bracket = 0, count_comma = 0;
    while(node < N_total){ 
      int startpos, endpos;
      std::string index_child1, index_child2, index_parent;
      std::string branch_length1, branch_length2;

      while(newick[i] == '(') i++;
      startpos = i;

      while(newick[i] != ':'){
        index_child1 += newick[i];
        i++;
      }
      i++;
      while(newick[i] != ',' && i < newick.size()){
        branch_length1 += newick[i];
        i++;
      }
      i++;
      if(newick[i] != '(' && i < newick.size()){
        while(newick[i] != ':'){
          index_child2 += newick[i];
          i++;
        } 
        i++;
        while(newick[i] != ')' && i < newick.size()){
          branch_length2 += newick[i];
          i++;
        }
        i++;
        endpos = i;

        //std::cerr << index_child1 << " " << index_child2 << " " << branch_length1 << " " << branch_length2 << std::endl;
        int child_left = stoi(index_child1), child_right = stoi(index_child2), parent = node;

        (*it_seq).tree.nodes[child_left].label  = child_left;
        (*it_seq).tree.nodes[child_right].label = child_right;
        (*it_seq).tree.nodes[parent].label      = parent;

        (*it_seq).tree.nodes[child_left].parent  = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[child_right].parent = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[child_left];
        (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[child_right];

        (*it_seq).tree.nodes[child_left].branch_length  = stof(branch_length1) * Ne;
        (*it_seq).tree.nodes[child_right].branch_length = stof(branch_length2) * Ne;

        //std::cerr << newick << std::endl;
        //std::cerr << newick.substr(startpos-1, endpos-startpos+1) << std::endl;
        newick.replace(startpos-1, endpos - startpos + 1, std::to_string(node)); 
        count_bracket = 0; 
        count_comma = 0;
        i = 0;
        for(; i < newick.size(); i++){
          if(newick[i] == '(') count_bracket++;
          if(newick[i] == ',') count_comma++;
        }
        //assert(count_comma == count_bracket);
        if(count_comma != count_bracket){
          break;
        }
        //std::cerr << newick << std::endl; 
        i = 0;
        node++;
      }
    }

    bool everyone_has_parent = true;
    for(int i = 0; i < N_total - 1; i++){
      if((*it_seq).tree.nodes[i].parent == NULL){
        everyone_has_parent = false;
        break;
      } 
    }
    //std::cerr << everyone_has_parent << std::endl;
    if(node != N_total || count_comma != count_bracket || !everyone_has_parent){
      seq.pop_back();
      num_tree--;
      it_seq--;
    }

    /*
       coordinates.resize(N_total); 
       (*it_seq).tree.GetCoordinates(coordinates);
       std::cerr << coordinates[N_total-1] << " " << coordinates[N_total-1]/Ne << std::endl; 
       */

    num_tree++;
    seq.emplace_back();
    it_seq = std::prev(seq.end(),1);
    //it_seq++;
  }

  seq.pop_back();

}


