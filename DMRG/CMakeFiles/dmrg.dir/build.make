# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_SOURCE_DIR = /home/yanir/Documents/cpp_stuff/DMRG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yanir/Documents/cpp_stuff/DMRG

# Include any dependencies generated for this target.
include CMakeFiles/dmrg.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/dmrg.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/dmrg.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dmrg.dir/flags.make

CMakeFiles/dmrg.dir/Mpo.cpp.o: CMakeFiles/dmrg.dir/flags.make
CMakeFiles/dmrg.dir/Mpo.cpp.o: Mpo.cpp
CMakeFiles/dmrg.dir/Mpo.cpp.o: CMakeFiles/dmrg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanir/Documents/cpp_stuff/DMRG/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dmrg.dir/Mpo.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dmrg.dir/Mpo.cpp.o -MF CMakeFiles/dmrg.dir/Mpo.cpp.o.d -o CMakeFiles/dmrg.dir/Mpo.cpp.o -c /home/yanir/Documents/cpp_stuff/DMRG/Mpo.cpp

CMakeFiles/dmrg.dir/Mpo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmrg.dir/Mpo.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanir/Documents/cpp_stuff/DMRG/Mpo.cpp > CMakeFiles/dmrg.dir/Mpo.cpp.i

CMakeFiles/dmrg.dir/Mpo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmrg.dir/Mpo.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanir/Documents/cpp_stuff/DMRG/Mpo.cpp -o CMakeFiles/dmrg.dir/Mpo.cpp.s

CMakeFiles/dmrg.dir/Main.cpp.o: CMakeFiles/dmrg.dir/flags.make
CMakeFiles/dmrg.dir/Main.cpp.o: Main.cpp
CMakeFiles/dmrg.dir/Main.cpp.o: CMakeFiles/dmrg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yanir/Documents/cpp_stuff/DMRG/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/dmrg.dir/Main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/dmrg.dir/Main.cpp.o -MF CMakeFiles/dmrg.dir/Main.cpp.o.d -o CMakeFiles/dmrg.dir/Main.cpp.o -c /home/yanir/Documents/cpp_stuff/DMRG/Main.cpp

CMakeFiles/dmrg.dir/Main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmrg.dir/Main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yanir/Documents/cpp_stuff/DMRG/Main.cpp > CMakeFiles/dmrg.dir/Main.cpp.i

CMakeFiles/dmrg.dir/Main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmrg.dir/Main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yanir/Documents/cpp_stuff/DMRG/Main.cpp -o CMakeFiles/dmrg.dir/Main.cpp.s

# Object files for target dmrg
dmrg_OBJECTS = \
"CMakeFiles/dmrg.dir/Mpo.cpp.o" \
"CMakeFiles/dmrg.dir/Main.cpp.o"

# External object files for target dmrg
dmrg_EXTERNAL_OBJECTS =

dmrg: CMakeFiles/dmrg.dir/Mpo.cpp.o
dmrg: CMakeFiles/dmrg.dir/Main.cpp.o
dmrg: CMakeFiles/dmrg.dir/build.make
dmrg: CMakeFiles/dmrg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yanir/Documents/cpp_stuff/DMRG/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable dmrg"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dmrg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dmrg.dir/build: dmrg
.PHONY : CMakeFiles/dmrg.dir/build

CMakeFiles/dmrg.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dmrg.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dmrg.dir/clean

CMakeFiles/dmrg.dir/depend:
	cd /home/yanir/Documents/cpp_stuff/DMRG && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yanir/Documents/cpp_stuff/DMRG /home/yanir/Documents/cpp_stuff/DMRG /home/yanir/Documents/cpp_stuff/DMRG /home/yanir/Documents/cpp_stuff/DMRG /home/yanir/Documents/cpp_stuff/DMRG/CMakeFiles/dmrg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dmrg.dir/depend

