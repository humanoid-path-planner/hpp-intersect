# Humanoid Path Planner - Intersect module

Copyright 2016 LAAS-CNRS

Author: Anna Seppala

## Description
HPP - INTERSECT is a library that implements intersection computation between fcl collision obstacles.

##Installation on ubuntu-14.04 64 bit with ros-indigo

To install HPP - INTERSECT, you will need to install one other package of the Humanoid Path Planner software with its respective dependecies. Please see the instructions below for the full installation of HPP - INTERSECT:

  1. Install HPP - FCL (make sure you are on branch "intersect" in the repository)
	- see https://github.com/anna-seppala/hpp-fcl

  2. Install Eigen 3
	- see http://eigen.tuxfamily.org/

  3. Clone the HPP - INTERSECT repository onto your local computer and update the submodule:

			git clone https://github.com/anna-seppala/hpp-intersect.git
			cd $HPP_INTERSECT_DIR/
			git submodule update --init --recursive

  4. Use CMake to install the HPP - INTERSECT library. For instance:

			mkdir build
			cd build
			cmake ..
			make install


##Documentation

Open $DEVEL_DIR/install/share/doc/hpp-intersect/doxygen-html/index.html in a web brower and you will have access to the code documentation.

