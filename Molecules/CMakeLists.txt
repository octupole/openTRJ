# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE AtomIndex.cpp  IteratorAtoms.cpp ResidueTypes.cpp
	Atoms.cpp  Groups.cpp  LCells.cpp  Topol.cpp  TopolPDB.cpp Percolation.cpp 
	Contacts.cpp Gyration.cpp GyrationJSON.cpp PercolationJSON.cpp FindCell.cpp)
add_library(mols STATIC ${SOURCE})
#target_include_directories(mols INTERFACE xdr/include voro++/include )
