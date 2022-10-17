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

// In this function, I will convert from anc/mut to tree sequence file format (functions are defined in ./include/src/tree_sequence.hpp)
void
ConvertToTreeSequence(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  bool compress = false;
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
  if(options.count("compress")){
    compress = true;
  }

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Converting " << options["anc"].as<std::string>() << " and " << options["mut"].as<std::string>() << " to tree sequence..." << std::endl;
  if (compress){
    std::cerr << "Combining equivalent branches and nodes; node age is stored as double precision in node metadata ..." << std::endl;
  }

  ////////// 1. Dump in tree sequence file format ////////
  
  std::string outfile = options["output"].as<std::string>() + ".trees";
  if (compress){
    //combine identical branches in adjacent trees, sticking average node time (weighted by span) into metadata
    DumpAsCompressedTreeSequence(options["anc"].as<std::string>(), options["mut"].as<std::string>(), outfile);
  } else {
    //new set of nodes for each tree (this does not compress trees)
    DumpAsTreeSequence(options["anc"].as<std::string>(), options["mut"].as<std::string>(), outfile);
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

// In this function, I will convert from anc/mut to tree sequence file format (functions are defined in ./include/src/tree_sequence.hpp)
void
ConvertFromTreeSequence(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("input")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, input." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Example code for converting from tree sequence file format." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Converting " << options["input"].as<std::string>() << " to " << options["anc"].as<std::string>() << " and " << options["mut"].as<std::string>() << "..." << std::endl;
  std::cerr << "Polytomies are broken at random (set seed)." << std::endl;

  bool no_branch_lengths = false;
  if(options.count("no_bl")){
    no_branch_lengths = true;
  }

  int seed = std::time(0) + getpid();
  if(options.count("seed")){
    seed = options["seed"].as<int>();
  } 

  ////////// 1. Dump in tree sequence file format ////////
  
  //new set of nodes for each tree (this does not compress trees)
  ConvertFromTreeSequence(options["anc"].as<std::string>(), options["mut"].as<std::string>(), options["input"].as<std::string>(), no_branch_lengths, seed);

  if(1){
  AncesTree anc;
  anc.Read(options["anc"].as<std::string>());
  anc.AssociateEquivalentBranches();
  anc.Dump(options["anc"].as<std::string>());
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

//////////////////////////////////////////////
int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("Conversion");
  options.add_options()
    ("help", "Print help.")
    ("compress", "Compress tree sequence by combining equivalent nodes in adjacent trees.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("anc", "Filename of file containing trees.", cxxopts::value<std::string>())
    ("mut", "Filename of file containing mut.", cxxopts::value<std::string>())
    ("no_bl", "If specified, assume that tree sequence has no branch lengths")
    ("seed", "Random seed (int).", cxxopts::value<int>())
    ("i,input", "Filename of input.", cxxopts::value<std::string>())
    ("o,output", "Filename of output (excl file extension).", cxxopts::value<std::string>());

  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("ConvertToTreeSequence")){

    ConvertToTreeSequence(options);

  }else if(!mode.compare("ConvertFromTreeSequence")){

    ConvertFromTreeSequence(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "ConvertToTreeSequence, ConvertFromTreeSequence." << std::endl;

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

