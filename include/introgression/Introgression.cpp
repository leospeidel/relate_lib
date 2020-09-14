#include "GetEquivalentBranches.cpp"

#include "cxxopts.hpp"
#include <string>

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("Relate");
  options.add_options()
    ("help", "Print help.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("i,input", "Filename of anc and mut files without file extension.", cxxopts::value<std::string>())
    ("o,output", "Filename for updated anc and mut files without file extension.", cxxopts::value<std::string>())
    ("ancient", "Filename of file containing bp of SNPs shared with ancient groups.", cxxopts::value<std::string>())
    ("branches", "Filename of file containing introgression branches (in .mut file format).", cxxopts::value<std::string>())
    ("poplabels", "Optional: Filename of file containing population labels. If ='hap', each haplotype is in its own group.", cxxopts::value<std::string>()) 
    ("num_bins", "Optional: Number of bins.", cxxopts::value<int>())
    ("mask", "Filename of file containing mask", cxxopts::value<std::string>())
    ("first_chr", "Optional: Index of fist chr", cxxopts::value<int>())
    ("last_chr", "Optional: Index of last chr", cxxopts::value<int>())
    ("seed", "Seed for MCMC in branch lengths estimation.", cxxopts::value<int>());

  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("GetEquivalentBranches")){

    GetEquivalentBranches(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "GetMutationsOnBranches." << std::endl;

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

