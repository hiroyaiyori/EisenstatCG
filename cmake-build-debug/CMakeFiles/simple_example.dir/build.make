# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/iyori/Desktop/UT/graduate/eisenstat

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/iyori/Desktop/UT/graduate/eisenstat/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/simple_example.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/simple_example.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simple_example.dir/flags.make

CMakeFiles/simple_example.dir/eisenstat.c.o: CMakeFiles/simple_example.dir/flags.make
CMakeFiles/simple_example.dir/eisenstat.c.o: ../eisenstat.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/iyori/Desktop/UT/graduate/eisenstat/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/simple_example.dir/eisenstat.c.o"
	/usr/local/Cellar/gcc/10.2.0/bin/gcc-10 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/simple_example.dir/eisenstat.c.o   -c /Users/iyori/Desktop/UT/graduate/eisenstat/eisenstat.c

CMakeFiles/simple_example.dir/eisenstat.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/simple_example.dir/eisenstat.c.i"
	/usr/local/Cellar/gcc/10.2.0/bin/gcc-10 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/iyori/Desktop/UT/graduate/eisenstat/eisenstat.c > CMakeFiles/simple_example.dir/eisenstat.c.i

CMakeFiles/simple_example.dir/eisenstat.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/simple_example.dir/eisenstat.c.s"
	/usr/local/Cellar/gcc/10.2.0/bin/gcc-10 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/iyori/Desktop/UT/graduate/eisenstat/eisenstat.c -o CMakeFiles/simple_example.dir/eisenstat.c.s

# Object files for target simple_example
simple_example_OBJECTS = \
"CMakeFiles/simple_example.dir/eisenstat.c.o"

# External object files for target simple_example
simple_example_EXTERNAL_OBJECTS =

simple_example: CMakeFiles/simple_example.dir/eisenstat.c.o
simple_example: CMakeFiles/simple_example.dir/build.make
simple_example: CMakeFiles/simple_example.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/iyori/Desktop/UT/graduate/eisenstat/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable simple_example"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simple_example.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simple_example.dir/build: simple_example

.PHONY : CMakeFiles/simple_example.dir/build

CMakeFiles/simple_example.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simple_example.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simple_example.dir/clean

CMakeFiles/simple_example.dir/depend:
	cd /Users/iyori/Desktop/UT/graduate/eisenstat/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/iyori/Desktop/UT/graduate/eisenstat /Users/iyori/Desktop/UT/graduate/eisenstat /Users/iyori/Desktop/UT/graduate/eisenstat/cmake-build-debug /Users/iyori/Desktop/UT/graduate/eisenstat/cmake-build-debug /Users/iyori/Desktop/UT/graduate/eisenstat/cmake-build-debug/CMakeFiles/simple_example.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/simple_example.dir/depend

