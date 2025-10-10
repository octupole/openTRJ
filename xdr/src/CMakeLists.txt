# Policy
cmake_minimum_required(VERSION 3.20) 

# Project name
#project(xdr-cmake C)

include_directories("${PROJECT_BINARY_DIR}/xdr")
add_library(xdrfile STATIC xdrfile_seek.c xdrfile_trr.c xdrfile_xtc.c xdrfile.c)
