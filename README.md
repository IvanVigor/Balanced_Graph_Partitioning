# Balanced Graph Partitioning
This is an implementation of the paper Balanced Graph Partitioning (authors: Konstantin Andreev, Harald Racke).  

## Libraries installation Windows 

The code was tested on python 3.5 with these additional libraries:
- numpy
- networkx
- python-igraph

In order to install them:

1) Open Python console 
2) write ` pip install numpy`  and wait termination
3) write ` pip install networkx`  and wait termination
4) write ` pip install python-igraph`  and wait termination
	In case of failure of this last step:
		1) go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#python-igraph
		2) write in python console: "pip install path/to/igraph.whl"

5) Then the file to launch is "Balancer_Cut"


		
## Library installation Linux 

The code wast tested on Python 2.7.11 because igraph library is not compatible with
Python 2.7+ versions.

### Requirements
On Linux (and other Unix-like systems), you will need a C and a C++ compiler, the tool 
"make" and the development header files for your version of Python. On Debian and Ubuntu 
Linux the "build-essential" and "python-dev" packages.

The code runs with these additional libraries:

- numpy
- networkx
- python-igraph

1) Open Python console 
2) write ` pip install numpy`  and wait termination
3) write ` pip install networkx`  and wait termination [remember to install the library in the right Python version (2.7)]
4) write ` pip install python-igraph`  and wait termination
	In case of failure of this last step:
		1) check "# Requirements"
		2) follow instructions on this page http://igraph.wikidot.com/installing-python-igraph-on-linux
		[Ubuntu Linux version tested by us]
		
5) Then the file to launch is "Balancer_Cut"
		
		
##  License
The code is provided with MIT license 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## References

Paper: "Balanced Graph Partitioning" http://www.math.cmu.edu/~kandreev/kpart.pdf

## Authors

[Ivan Vigorito](https://github.com/IvanVigor) and Lorenzo Frigerio
