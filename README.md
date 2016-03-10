# elss

**elss** is a Monte Carlo (MC) code for the analysis and the inference of energy landscapes in protein sequence spaces. 
**elss** can either be used to simulate a trajectory with a user-defined energy function, or to infer a data-driven statistical model from a multiple sequence alignment (MSA) 
via maximum a posteriori (MAP) estimation. **elss** implements an energy function given by a sum of potential terms that control both the frequencies of amino acids 
at the different positions along the chain and their correlations. 

# A simple example

    $ elss -p prm -n 10000 

the program will read custom values for the parameters of the energy function from file _prm_ (**-p prm**); then it will simulate a (10000 sweeps long) trajectory (**-n 10000**) starting from a random sequence. Sampled sequences will be saved to file **0.trj** (in FASTA format) every 10 sweeps. 

    $ cat 0.trj
    >       0     -113.204       -4.367     -108.837
    IYVGNLPYTSTHEDLNTHAKTYDEIENVHMAYSD-SNFRGFAFVEFHDKEDAAKALSG--D---------
    >      10     -121.961       -4.167     -117.794
    IYIGNLIYSMTNDELTQAFETYGDISRVSIILDDTGKPKGYAFVRFVEKEGVKLCVAQLKA---------
    >      20     -100.255       -3.609      -96.646
    LYMDLLDESVSAADLKVQFGKFGEECRVFHVRDANGFSKGRAFVEFKSRDQAVSAMEHR--R--------
    ...
    >   10000      -96.654       -4.599      -92.055
    VFVGGLPQSVTTA-LLEHFTDAGEKEEVMIGGDDSERSKGFAFVTFCQEEQCEKFVDESNYKEILGRQV-
    
Same analysis with small differences: 

    $ elss -p prm --seq start.fa -n 10000 --nupdate 1 
    
In this case, the starting sequence will be read from file **start.fa** (**--seq start.fa**) in FASTA format, and sequences will be dumped at every sweep (**--nupdate 1**). 

# A slightly more complex example

    $ mpiexec -n 8 elss --fasta PF00076.fa --learn-agd 2000 -n 10000 --lambda 100.0
    $ elss -r rst -n 100000

First line: the program reads a MSA from file **PF00076.fa** (**--fasta PF00076.fa**) and compute the maximum-a-posteriori estimate of the parameters of the energy function. The algorithm takes 2000 accelerated gradient descent steps (**--learn-agd 2000**), computing the gradient of the objective function from the accumulated statistics of 8 (**mpiexec -n 8**) _persistent_, 10000 sweeps long Markov chains (**--nsweeps 10000**). 
The option --lambda controls the value of the hyper-parameter of the Gaussian prior, P_0(\theta) $\propto$ \exp (-0.5 * \lambda * \theta^2). Higher values correspond to more strongly regularized solutions. The default is 100.0. 
Output files: the files _prm_ (a file containing the estimated parameters), _rst_ (a binary restart file) and _LEARN.log_, a log file. 

Second line: **elss** will read the final energy model from the previous calculation (_-r rst_) and simulate a (100000 sweeps long) trajectory in sequence space starting from a random sequence. The command will dump the file _0.trj__ (see previous example). 

# Format of the parameter file
The parameter file is the result of the inference algorithm, and can be modified to investigate the effect of perturbations to the energy function. 

    $ cat prm
    # flag L q
    protein 70 21
    # biases
    1  2 -2.0 
    20 21  0.5
    ...
    # couplings
    55 70  5  20 1.0
    ...

- the \# symbol is used to mark comments
- the first line contains: 
    - a flag describing the output format ("_protein_", do not change this)
    - the number of positions in the protein and 
    - the number of symbols in the amino acids alphabet (21 = 20 natural amino acids + the gap symbol (\-)
- the next lines define the parameters of the energy function
- amino acids are coded using their numerical index in alphabetical order (gap has index 21)
- lines with three fields define specific biases e.g.:  
    the line "1 1 -2.0" introduces a potential term in the energy function that stabilizes (-2.0)
    the amino acid no. 2 (CYS) at position 1. 
- lines with five fields define interaction terms between amino acids at different positions.  
    the line "55 70 5 20  1.0" in the example introduces a potential term in the energy function that disfavours    
    the amino acid pair PHE and TYR at positions no. 55 and 70 respectevly. 

# Prerequisites

In order to compile **elss**, you will need to have a **Fortran compiler** installed on your machine.   
For gfortran, it is necessary to use version 4.6.0 or above.
If you are using Debian or a Debian derivative such as Ubuntu, you can install the gfortran compiler using the following command:

    sudo apt-get install gfortran

The inference algorithm works by simulating a swarm of persistent Markov chains. 
To compile **elss** with support for parallel runs on a distributed-memory architecture, you will need to have a valid **MPI implementation** installed on your machine. 
The code has been tested and is known to work with the latest versions of both OpenMPI and MPICH. 
OpenMPI (recommended) can be installed on Debian derivatives with:

    sudo apt-get install openmpi-bin libopenmpi1.6 libopenmpi-dev

For details on running MPI jobs with OpenMPI see [this link](https://www.open-mpi.org/faq/?category=running)

Alternatively, MPICH can be installed like this:

    sudo apt-get install mpich libmpich-dev

The compiling and linking of source files is handled by **Gnu Make**. 
If you are using Debian or a Debian derivative such as Ubuntu, you should find 3.81 already installed. 

(optional) **git** version control software for obtaining source code

    sudo apt-get install git

# Obtaining the source

All elss source code is hosted on GitHub. You can download the source code directly from GitHub or, if you have 
the git version control software installed on your computer, you can use git to obtain the source code. 
The latter method has the benefit that it is easy to receive updates directly from the GitHub repository. 
GitHub has a good set of instructions for how to set up git to work with GitHub since this involves setting up ssh keys. 
With git installed and setup, the following command will download the full source code from the GitHub repository:

    git clone https://github.com/yyy/xxx.git

# Compiling on Linux and Mac OS X ###

To compile elss on Linux or Max OS X, from the project root directory enter the src directory and type make:

    (cd src; make) 

This will build the executable elss.

# Testing

Run the _run-test.bash_ script in the test directory: 

    (cd test; ./run-test.bash)

# Contributing

**elss** is an OPEN Source Project so please help out by [reporting bugs](http://github.com/simomarsili/elss/issues) or [forking and opening pull](https://github.com/simomarsili/elss) requests when possible.

# LICENSE (MIT/X11)

Copyright (c) 2016 Simone Marsili ("Author")

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

