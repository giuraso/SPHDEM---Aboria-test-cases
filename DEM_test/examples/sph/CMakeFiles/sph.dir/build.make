# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/giuseppe/code/Aboria/examples/sph

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/giuseppe/code/Aboria/examples/sph

# Include any dependencies generated for this target.
include CMakeFiles/sph.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/sph.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/sph.dir/flags.make

CMakeFiles/sph.dir/sph.cpp.o: CMakeFiles/sph.dir/flags.make
CMakeFiles/sph.dir/sph.cpp.o: sph.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/giuseppe/code/Aboria/examples/sph/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/sph.dir/sph.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/sph.dir/sph.cpp.o -c /home/giuseppe/code/Aboria/examples/sph/sph.cpp

CMakeFiles/sph.dir/sph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sph.dir/sph.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/giuseppe/code/Aboria/examples/sph/sph.cpp > CMakeFiles/sph.dir/sph.cpp.i

CMakeFiles/sph.dir/sph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sph.dir/sph.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/giuseppe/code/Aboria/examples/sph/sph.cpp -o CMakeFiles/sph.dir/sph.cpp.s

CMakeFiles/sph.dir/sph.cpp.o.requires:
.PHONY : CMakeFiles/sph.dir/sph.cpp.o.requires

CMakeFiles/sph.dir/sph.cpp.o.provides: CMakeFiles/sph.dir/sph.cpp.o.requires
	$(MAKE) -f CMakeFiles/sph.dir/build.make CMakeFiles/sph.dir/sph.cpp.o.provides.build
.PHONY : CMakeFiles/sph.dir/sph.cpp.o.provides

CMakeFiles/sph.dir/sph.cpp.o.provides.build: CMakeFiles/sph.dir/sph.cpp.o

# Object files for target sph
sph_OBJECTS = \
"CMakeFiles/sph.dir/sph.cpp.o"

# External object files for target sph
sph_EXTERNAL_OBJECTS =

sph: CMakeFiles/sph.dir/sph.cpp.o
sph: /usr/lib/libvtkCommon.so.5.8.0
sph: /usr/lib/libvtkFiltering.so.5.8.0
sph: /usr/lib/libvtkImaging.so.5.8.0
sph: /usr/lib/libvtkGraphics.so.5.8.0
sph: /usr/lib/libvtkGenericFiltering.so.5.8.0
sph: /usr/lib/libvtkIO.so.5.8.0
sph: /usr/lib/libvtkRendering.so.5.8.0
sph: /usr/lib/libvtkVolumeRendering.so.5.8.0
sph: /usr/lib/libvtkHybrid.so.5.8.0
sph: /usr/lib/libvtkWidgets.so.5.8.0
sph: /usr/lib/libvtkParallel.so.5.8.0
sph: /usr/lib/libvtkInfovis.so.5.8.0
sph: /usr/lib/libvtkGeovis.so.5.8.0
sph: /usr/lib/libvtkViews.so.5.8.0
sph: /usr/lib/libvtkCharts.so.5.8.0
sph: /usr/lib/libvtkViews.so.5.8.0
sph: /usr/lib/libvtkInfovis.so.5.8.0
sph: /usr/lib/libvtkWidgets.so.5.8.0
sph: /usr/lib/libvtkVolumeRendering.so.5.8.0
sph: /usr/lib/libvtkHybrid.so.5.8.0
sph: /usr/lib/libvtkParallel.so.5.8.0
sph: /usr/lib/libvtkRendering.so.5.8.0
sph: /usr/lib/libvtkImaging.so.5.8.0
sph: /usr/lib/libvtkGraphics.so.5.8.0
sph: /usr/lib/libvtkIO.so.5.8.0
sph: /usr/lib/libvtkFiltering.so.5.8.0
sph: /usr/lib/libvtkCommon.so.5.8.0
sph: /usr/lib/libvtksys.so.5.8.0
sph: CMakeFiles/sph.dir/build.make
sph: CMakeFiles/sph.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable sph"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sph.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/sph.dir/build: sph
.PHONY : CMakeFiles/sph.dir/build

CMakeFiles/sph.dir/requires: CMakeFiles/sph.dir/sph.cpp.o.requires
.PHONY : CMakeFiles/sph.dir/requires

CMakeFiles/sph.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sph.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sph.dir/clean

CMakeFiles/sph.dir/depend:
	cd /home/giuseppe/code/Aboria/examples/sph && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/giuseppe/code/Aboria/examples/sph /home/giuseppe/code/Aboria/examples/sph /home/giuseppe/code/Aboria/examples/sph /home/giuseppe/code/Aboria/examples/sph /home/giuseppe/code/Aboria/examples/sph/CMakeFiles/sph.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/sph.dir/depend
