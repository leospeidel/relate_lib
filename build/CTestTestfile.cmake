# CMake generated Testfile for 
# Source directory: /data/smew1/speidel/genomics/relate_lib
# Build directory: /data/smew1/speidel/genomics/relate_lib/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest "/data/smew1/speidel/genomics/relate_lib/bin/Tests")
set_tests_properties(UnitTest PROPERTIES  _BACKTRACE_TRIPLES "/data/smew1/speidel/genomics/relate_lib/CMakeLists.txt;43;add_test;/data/smew1/speidel/genomics/relate_lib/CMakeLists.txt;0;")
subdirs("include/src")
subdirs("include/src/gzstream")
subdirs("include/src/tskit")
subdirs("include/test")
subdirs("include/example")
