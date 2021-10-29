#include "catch_amalgamated.hpp"

#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "mutations.hpp"
#include "sample.hpp"


TEST_CASE( "Test parsing anc/mut files, trees" ){

  std::string filename_anc = "../include/test/data/output_files/example.anc";
  std::string filename_mut = "../include/test/data/output_files/example.mut";

  AncMutIterators ancmut(filename_anc, filename_mut);

  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0;

  //iterate through whole file
  int count_trees = 0;
  while(ancmut.NextTree(mtr, it_mut) >= 0.0){
    count_trees++;
  }
  REQUIRE((*it_mut).pos == 3999486); //it_mut points to where last tree starts
  REQUIRE(count_trees == 1366); //total tree count should match;
  //REQUIRE(ancmut.NextSNP(mtr, it_mut) == -1.0);

  ancmut.OpenFiles(filename_anc, filename_mut);
  ancmut.CloseFiles();
  ancmut.CloseFiles();

}

TEST_CASE( "Test parsing anc/mut files, SNPs" ){

  std::string filename_anc = "../include/test/data/output_files/example.anc";
  std::string filename_mut = "../include/test/data/output_files/example.mut";

  AncMutIterators ancmut(filename_anc, filename_mut);

  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0;

  //get first SNP (Only necessary when using NextSNP, for NextTree I don't need to call this)
  ancmut.FirstSNP(mtr, it_mut);
  REQUIRE((*it_mut).pos == 73);
  REQUIRE(mtr.pos == 0);

  //iterate through whole file
  int count_SNP = 1; //first tree already loaded
  while(ancmut.NextSNP(mtr, it_mut) >= 0.0){
    count_SNP++;
  }
  REQUIRE(count_SNP == ancmut.NumSnps()); //total tree count should match;
  
  ancmut.CloseFiles();

}

