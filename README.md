# atINFLATE_2.0
aaRS-tRNA Interaction Network Fitness LAndscape Express 2.0 with rate selection parameter

# atINFLATE

**atINFLATE: aaRS-tRNA Interaction Network Fitness Landscape
Topographer Express**

**Author:** Andrea Collins-Hed

atINFLATE simulates a model for the origin of the aaRS-tRNA
interaction network in coevolution with a large set of protein-coding
genes to express a target proteome. 

## Getting Started

Clone the atINFLATE git repository and change directory into your local
clone to compile atINFLATE

## Pre-Installation of Dependencies

Compiling atINFLATE requires a local copy of version 3.3.9 of the
[Eigen C++ template library for linear
algebra](https://eigen.tuxfamily.org/dox/GettingStarted.html). 

atINFLATE uses the KroneckerProduct module from the `unsupported`
collection of modules of Eigen.

## Compiling

atINFLATE uses one g++ specific instruction (the
"__builtin_popcount()" instruction which outputs the weight of a
number in the coding sense) so it must be complied with g++. 

The command to compile is:

`g++ -I /path/to/eigen/ atinflate.cpp -std=c++11 -O2 -o atinflate`

Or, on Linux and Mac OS X, if you symlink or copy the Eigen folder into
/usr/local/include/, you can compile with: 

`g++ atinflate.cpp -std=c++11 -O2 -o atinflate`

Although currently atINFLATE does not yet utilize any openMP code, openMP may
become a dependency in the future, in which case the `-fopenmp` option
should be added to the command to compile. 

## Help with Running 

Once compiled, `atinflate --help` will output a help page similar to a
`man` page.
