# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /data/smew1/speidel/genomics/relate_lib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /data/smew1/speidel/genomics/relate_lib/build

# Include any dependencies generated for this target.
include include/src/gzstream/CMakeFiles/gzstreamShared.dir/depend.make

# Include the progress variables for this target.
include include/src/gzstream/CMakeFiles/gzstreamShared.dir/progress.make

# Include the compile flags for this target's objects.
include include/src/gzstream/CMakeFiles/gzstreamShared.dir/flags.make

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o: include/src/gzstream/CMakeFiles/gzstreamShared.dir/flags.make
include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o: ../include/src/gzstream/gzstream.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/smew1/speidel/genomics/relate_lib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o"
	cd /data/smew1/speidel/genomics/relate_lib/build/include/src/gzstream && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gzstreamShared.dir/gzstream.cpp.o -c /data/smew1/speidel/genomics/relate_lib/include/src/gzstream/gzstream.cpp

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gzstreamShared.dir/gzstream.cpp.i"
	cd /data/smew1/speidel/genomics/relate_lib/build/include/src/gzstream && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/smew1/speidel/genomics/relate_lib/include/src/gzstream/gzstream.cpp > CMakeFiles/gzstreamShared.dir/gzstream.cpp.i

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gzstreamShared.dir/gzstream.cpp.s"
	cd /data/smew1/speidel/genomics/relate_lib/build/include/src/gzstream && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/smew1/speidel/genomics/relate_lib/include/src/gzstream/gzstream.cpp -o CMakeFiles/gzstreamShared.dir/gzstream.cpp.s

# Object files for target gzstreamShared
gzstreamShared_OBJECTS = \
"CMakeFiles/gzstreamShared.dir/gzstream.cpp.o"

# External object files for target gzstreamShared
gzstreamShared_EXTERNAL_OBJECTS =

../bin/libgzstreamShared.so: include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o
../bin/libgzstreamShared.so: include/src/gzstream/CMakeFiles/gzstreamShared.dir/build.make
../bin/libgzstreamShared.so: /usr/lib64/libz.so
../bin/libgzstreamShared.so: include/src/gzstream/CMakeFiles/gzstreamShared.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/data/smew1/speidel/genomics/relate_lib/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../../bin/libgzstreamShared.so"
	cd /data/smew1/speidel/genomics/relate_lib/build/include/src/gzstream && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gzstreamShared.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
include/src/gzstream/CMakeFiles/gzstreamShared.dir/build: ../bin/libgzstreamShared.so

.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/build

include/src/gzstream/CMakeFiles/gzstreamShared.dir/clean:
	cd /data/smew1/speidel/genomics/relate_lib/build/include/src/gzstream && $(CMAKE_COMMAND) -P CMakeFiles/gzstreamShared.dir/cmake_clean.cmake
.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/clean

include/src/gzstream/CMakeFiles/gzstreamShared.dir/depend:
	cd /data/smew1/speidel/genomics/relate_lib/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/smew1/speidel/genomics/relate_lib /data/smew1/speidel/genomics/relate_lib/include/src/gzstream /data/smew1/speidel/genomics/relate_lib/build /data/smew1/speidel/genomics/relate_lib/build/include/src/gzstream /data/smew1/speidel/genomics/relate_lib/build/include/src/gzstream/CMakeFiles/gzstreamShared.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/depend

