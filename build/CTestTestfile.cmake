# CMake generated Testfile for 
# Source directory: /Users/leo/Documents/genomics/relate_lib
# Build directory: /Users/leo/Documents/genomics/relate_lib/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest "/Users/leo/Documents/genomics/relate_lib/bin/Tests")
subdirs(include/src)
subdirs(include/src/gzstream)
subdirs(include/src/tskit)
subdirs(include/test)
subdirs(include/example)
