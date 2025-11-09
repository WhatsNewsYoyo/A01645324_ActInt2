# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( Quinto )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )

find_package(CGAL REQUIRED COMPONENTS Qt5)
add_definitions(-DCGAL_USE_BASIC_VIEWER)

create_single_source_cgal_program( "cgal.cpp" )
create_single_source_cgal_program( "voronoi.cpp" )
