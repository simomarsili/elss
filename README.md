# mcsg

**mcsg** is a Monte Carlo (MC) code for the analysis and the inference of energy landscapes in protein sequence spaces,
within the framework of [Direct Coupling Analysis (DCA)](https://en.wikipedia.org/wiki/Direct_coupling_analysis).
**mcsg** can either be used to simulate a trajectory with a user-defined energy function, or to infer a data-driven statistical model
from a multiple sequence alignment (MSA) via maximum a posteriori (MAP) estimation. Learning of the model's parameters is carried out using a
stochastic approximation of the gradients of the partition function based on Markov chain Monte Carlo simulations. 

# Obtaining the source

All **mcsg** source code is hosted on Github. 
You can download the latest version of the code using [this link](https://github.com/simomarsili/mcsg/archive/v0.3.1.tar.gz). 

# Prerequisites

In order to compile **mcsg**, you will need to have a **Fortran compiler** installed on your machine.   
For gfortran, it is necessary to use version 4.6.0 or above.
If you are using Debian or a Debian derivative such as Ubuntu, you can install the gfortran compiler using the following command:

    sudo apt-get install gfortran

The inference algorithm works by simulating a swarm of persistent Markov chains. 
To compile **mcsg** with support for parallel runs on a distributed-memory architecture, you will need to have a valid **MPI implementation** installed on your machine. 
The code has been tested and is known to work with the latest versions of both OpenMPI and MPICH.   
OpenMPI (recommended) can be installed on Debian derivatives with:

    sudo apt-get install openmpi-bin libopenmpi1.10 libopenmpi-dev

For details on running MPI jobs with OpenMPI see [this link](https://www.open-mpi.org/faq/?category=running)

Alternatively, MPICH can be installed with:

    sudo apt-get install mpich libmpich-dev

The compiling and linking of source files is handled by **Gnu Make**. 
If you are using Debian or a Debian derivative such as Ubuntu, you should find 3.81 already installed. 

(optional) **git** version control software for obtaining source code

    sudo apt-get install git

# Compiling ###

To compile **mcsg**, from the project root directory enter the src directory and type make:

    (cd src; make) 

This will build the executable mcsg.

# Testing

Run the **run-test.bash** script in the test directory: 

    (cd test; ./run-test.bash)

# A simple example

    $ mcsg -p prm -n 10000 

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

    $ mcsg -p prm --seq start.fa -n 10000 --nupdate 1 
    
In this case, the starting sequence will be read from file **start.fa** (**--seq start.fa**) in FASTA format, and sequences will be dumped at every sweep (**--nupdate 1**). 

# A slightly more complex example

    $ mpiexec -n 8 mcsg --fasta PF00076.fa --learn-agd 2000 -n 10000 --lambda 0.01
    $ mcsg -r rst -n 100000

First line: the program reads a MSA from file **PF00076.fa** (**--fasta PF00076.fa**) and compute the maximum-a-posteriori estimate of the parameters of the energy function. The algorithm takes 2000 accelerated gradient descent steps (**--learn-agd 2000**), computing the gradient of the objective function from the accumulated statistics of 8 (**mpiexec -n 8**) _persistent_, 10000 sweeps long Markov chains (**--nsweeps 10000**). 
The option --lambda controls regularization strength. Higher values correspond to more regularized solutions. The default is 0.01. 
Output files: the files _prm_ (a file containing the estimated parameters), _rst_ (a binary restart file) and _LEARN.log_, a log file. 

Second line: **mcsg** will read the final energy model from the previous calculation (_-r rst_) and simulate a (100000 sweeps long) trajectory in sequence space starting from a random sequence. The command will dump the file _0.trj__ (see previous example). 

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

# Contributing

**mcsg** is an OPEN Source Project so please help out by [reporting bugs](http://github.com/simomarsili/mcsg/issues) or [forking and opening pull](https://github.com/simomarsili/mcsg) requests when possible.

# LICENSE (BSD 3 clause)

Copyright (c) 2016, Simone Marsili 
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
