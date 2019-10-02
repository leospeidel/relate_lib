#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "mutations.hpp"
#include "sample.hpp"
#include "tree_sequence.hpp"

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

// In this function, I list 3 different ways of parsing anc/mut files
void
Parse(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Example code for parsing anc/mut files." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Parsing " << options["anc"].as<std::string>() << " and " << options["mut"].as<std::string>() << "..." << std::endl;

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;


  ////////// 1. Read one tree at a time /////////

  //We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
  //The mut file is read once, file is closed after constructor is called.
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());

  //iterate through whole file
  while(num_bases_tree_persists >= 0.0){

    //mtr stores the marginal tree, it_mut points to the first SNP at which the marginal tree starts.
    //If marginal tree has no SNPs mapping to it, it_mut points to the first SNP of the following tree that has a SNP mapping to it
    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

    //mut can be modified, it_mut can be incremented etc. 
    //However, mut is used in ancmut, so any changes may affect how trees and snps are parsed.
    (*it_mut).pos = 1000;

  }
  //at this point, anc file has been closed.


  ////////// 2. Iterate through all SNPs /////////

  float num_bases_SNP_persists = 0.0;

  //We can now reopen the files (if anc file was still open, this will automatically close the file first)
  ancmut.OpenFiles(options["anc"].as<std::string>(), options["mut"].as<std::string>());

  //Let's now read one SNP at a time:

  //get first SNP (Only necessary when using NextSNP, for NextTree I don't need to call this)
  num_bases_SNP_persists = ancmut.FirstSNP(mtr, it_mut);

  //iterate through whole file
  while(num_bases_SNP_persists >= 0.0){

    //it_mut points to a SNP, mtr stores the marginal tree corresponding to that SNP.
    num_bases_SNP_persists = ancmut.NextTree(mtr, it_mut);

    //mut can be modified, it_mut can be incremented etc. 
    //However, mut is used in ancmut, so any changes may affect how trees and snps are parsed.
    (*it_mut).pos = 1000;

  }
  //at this point, anc file has been closed.


  ////////// 3. Simply read anc/mut all at once /////////

  //For small examples, we can read the entire anc file at once.

  AncesTree anc; //anc.seq is a list of MarginalTrees (type CorrTrees)
  Mutations mut; //mut.info is a vector of SNPInfo (type Muts)

  anc.Read(options["anc"].as<std::string>());
  mut.Read(options["mut"].as<std::string>());

  //anc for dumping these files:
  anc.Dump(options["output"].as<std::string>() + ".anc");
  mut.Dump(options["output"].as<std::string>() + ".mut");

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

//In this function, I am parsing haps/sample and poplabels files
void
ParseData(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("haps") || !options.count("sample") || !options.count("poplabels")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, sample, poplabels." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Example code for parsing haps/sample and poplabels files." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Parsing " << options["haps"].as<std::string>() << ", " << options["sample"].as<std::string>() << ", and " << options["poplabels"].as<std::string>() << "..." << std::endl;


  haps m_hap(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  Data data(m_hap.GetN(), m_hap.GetL());

  std::vector<char> sequence(data.N);
  int bp;
  //parse snp by snp
  for(int snp = 0; snp < data.L; snp++){ 
    m_hap.ReadSNP(sequence, bp);
  }

  //Parse poplabels file
  //(defined in include/src/sample.hpp)
  Sample poplabels;
  poplabels.Read(options["poplabels"].as<std::string>());
  poplabels.AssignPopOfInterest("group1,group2"); //assign populations of interest (comma separated)

  //group of haplotype 0:
  std::string group  = poplabels.groups[poplabels.group_of_haplotype[0]];
  //region of group
  std::string region = poplabels.region_of_group[group];

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

// In this function, I will convert from anc/mut to tree sequence file format (functions are defined in ./include/src/tree_sequence.hpp)
void
ConvertToTreeSequence(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Example code for converting to tree sequence file format." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Converting " << options["anc"].as<std::string>() << " and " << options["mut"].as<std::string>() << "  to tree sequence..." << std::endl;

  ////////// 1. Dump in tree sequence file format ////////
  
  //new set of nodes for each tree (this does not compress trees)
  DumpAsTreeSequence(options["anc"].as<std::string>(), options["mut"].as<std::string>(), options["output"].as<std::string>() + ".trees");

  //ignore branch lengths and combine identical branches in adjacent trees
  DumpAsTreeSequenceTopoOnly(options["anc"].as<std::string>(), options["mut"].as<std::string>(), options["output"].as<std::string>() + ".topo.trees");

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

//////////////////////////////////////////////
int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("RelateExtract");
  options.add_options()
    ("help", "Print help.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("poplabels", "Filename of file containing population labels.", cxxopts::value<std::string>()) 
    ("anc", "Filename of file containing trees.", cxxopts::value<std::string>())
    ("mut", "Filename of file containing mut.", cxxopts::value<std::string>())
    ("haps", "Filename of haps file (Output file format of Shapeit).", cxxopts::value<std::string>())
    ("sample", "Filename of sample file (Output file format of Shapeit).", cxxopts::value<std::string>())
    ("pop_of_interest", "Population label. If not specified, use all haplotypes.", cxxopts::value<std::string>())
    ("bp_of_interest", "BP of position of interest.", cxxopts::value<int>())
    ("first_bp", "BP of first SNP of interest.", cxxopts::value<int>())
    ("last_bp", "BP of last SNP of interest.", cxxopts::value<int>())
    ("o,output", "Filename of output (excl file extension).", cxxopts::value<std::string>());

  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("Parse")){

    Parse(options);

  }else if(!mode.compare("ParseData")){

    ParseData(options);

  }else if(!mode.compare("ConvertToTreeSequence")){

    ConvertToTreeSequence(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "Parse, ParseData, ConvertToTreeSequence." << std::endl;

  }

  bool help = false;
  if(!options.count("mode")){
    std::cout << "Not enough arguments supplied." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    exit(0);
  }

  return 0;

}

