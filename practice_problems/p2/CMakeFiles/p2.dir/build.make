# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kvasudev/FEM_practice

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kvasudev/FEM_practice

# Include any dependencies generated for this target.
include practice_problems/p2/CMakeFiles/p2.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include practice_problems/p2/CMakeFiles/p2.dir/compiler_depend.make

# Include the progress variables for this target.
include practice_problems/p2/CMakeFiles/p2.dir/progress.make

# Include the compile flags for this target's objects.
include practice_problems/p2/CMakeFiles/p2.dir/flags.make

practice_problems/p2/CMakeFiles/p2.dir/p2.cc.o: practice_problems/p2/CMakeFiles/p2.dir/flags.make
practice_problems/p2/CMakeFiles/p2.dir/p2.cc.o: practice_problems/p2/p2.cc
practice_problems/p2/CMakeFiles/p2.dir/p2.cc.o: practice_problems/p2/CMakeFiles/p2.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/kvasudev/FEM_practice/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object practice_problems/p2/CMakeFiles/p2.dir/p2.cc.o"
	cd /home/kvasudev/FEM_practice/practice_problems/p2 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT practice_problems/p2/CMakeFiles/p2.dir/p2.cc.o -MF CMakeFiles/p2.dir/p2.cc.o.d -o CMakeFiles/p2.dir/p2.cc.o -c /home/kvasudev/FEM_practice/practice_problems/p2/p2.cc

practice_problems/p2/CMakeFiles/p2.dir/p2.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/p2.dir/p2.cc.i"
	cd /home/kvasudev/FEM_practice/practice_problems/p2 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kvasudev/FEM_practice/practice_problems/p2/p2.cc > CMakeFiles/p2.dir/p2.cc.i

practice_problems/p2/CMakeFiles/p2.dir/p2.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/p2.dir/p2.cc.s"
	cd /home/kvasudev/FEM_practice/practice_problems/p2 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kvasudev/FEM_practice/practice_problems/p2/p2.cc -o CMakeFiles/p2.dir/p2.cc.s

# Object files for target p2
p2_OBJECTS = \
"CMakeFiles/p2.dir/p2.cc.o"

# External object files for target p2
p2_EXTERNAL_OBJECTS =

practice_problems/p2/p2: practice_problems/p2/CMakeFiles/p2.dir/p2.cc.o
practice_problems/p2/p2: practice_problems/p2/CMakeFiles/p2.dir/build.make
practice_problems/p2/p2: /usr/local/lib/libdeal_II.so.9.7.0-pre
practice_problems/p2/p2: /usr/lib64/libtbb.so
practice_problems/p2/p2: /usr/lib64/libz.so
practice_problems/p2/p2: /usr/lib64/libumfpack.so
practice_problems/p2/p2: /usr/lib64/libcholmod.so
practice_problems/p2/p2: /usr/lib64/libccolamd.so
practice_problems/p2/p2: /usr/lib64/libcolamd.so
practice_problems/p2/p2: /usr/lib64/libcamd.so
practice_problems/p2/p2: /usr/lib64/libsuitesparseconfig.so
practice_problems/p2/p2: /usr/lib64/libamd.so
practice_problems/p2/p2: /usr/lib64/libopenblas.so
practice_problems/p2/p2: /usr/lib64/libmetis.so
practice_problems/p2/p2: /usr/lib64/libmuparser.so
practice_problems/p2/p2: practice_problems/p2/CMakeFiles/p2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/kvasudev/FEM_practice/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable p2"
	cd /home/kvasudev/FEM_practice/practice_problems/p2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/p2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
practice_problems/p2/CMakeFiles/p2.dir/build: practice_problems/p2/p2
.PHONY : practice_problems/p2/CMakeFiles/p2.dir/build

practice_problems/p2/CMakeFiles/p2.dir/clean:
	cd /home/kvasudev/FEM_practice/practice_problems/p2 && $(CMAKE_COMMAND) -P CMakeFiles/p2.dir/cmake_clean.cmake
.PHONY : practice_problems/p2/CMakeFiles/p2.dir/clean

practice_problems/p2/CMakeFiles/p2.dir/depend:
	cd /home/kvasudev/FEM_practice && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kvasudev/FEM_practice /home/kvasudev/FEM_practice/practice_problems/p2 /home/kvasudev/FEM_practice /home/kvasudev/FEM_practice/practice_problems/p2 /home/kvasudev/FEM_practice/practice_problems/p2/CMakeFiles/p2.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : practice_problems/p2/CMakeFiles/p2.dir/depend

