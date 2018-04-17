# openTRJ

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a2d6bcff43164c28b7f8c147334b4065)](https://app.codacy.com/app/octupole/openTRJ?utm_source=github.com&utm_medium=referral&utm_content=octupole/openTRJ&utm_campaign=badger)

## Synopsis

openTRJ is an ensemble of programs to compute a few properties from MD simulation trajectories written in .dcd or .xtc format. It was created 
to treat systems undergoing self-aggregations, such as micelles, inverted micelles, aggregates of proteins. The compilation with no parameters 
create three programs:
* trjVoronoi which computes the Voronoi volumes using voro++
* trjProp    which computes radius of gyration and aggregation properties
* trjSaxs    which computes the Saxs profiles of the system

All programs require a .pdb file containing the coordinates of the whole system, including waters and ions if needed. trjVoronoi uses voro++ 
by Chris H. Rycroft to compute properties of molecular systems in solution from molecular dynamics trajectories 
obtained with GROMACS and NAMD. 

## Motivation

This project was created to compute the self aggregation properties of large systems.

## Installation

openTRJ uses standard cmake (version > 3.0) scripts to generate Makefiles. This means that a typical installation will consist in creating a build directory and running within
that directory a command such as 'cmake ..' followed by <br/>

make<br />
make install<br />

The code needs that the fftw3 libraries be installed on the system. If openMPI and opneMP are available the code is compiled for parallel runs.

openTRJ uses several features available with the C++11 standard, thus it will not compile with GNU g++ compilers earlier 
than 5.4. It has been successfully compiled on Centos, Debian and 
Mac OS X systems supporting those versions of the C++ compilers.

## How to cite openTRJ: 

To be added

## License

  CeCILL FREE SOFTWARE LICENSE AGREEMENT

Version 2.1 dated 2013-06-21


This Agreement is a Free Software license agreement that is the result
of discussions between its authors in order to ensure compliance with
the two main principles guiding its drafting:

  * firstly, compliance with the principles governing the distribution
    of Free Software: access to source code, broad rights granted to users,
  * secondly, the election of a governing law, French law, with which it
    is conformant, both as regards the law of torts and intellectual
    property law, and the protection that it offers to both authors and
    holders of the economic rights over software.

Find the full licence here: http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.txt
