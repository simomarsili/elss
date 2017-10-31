# **elss**

**elss** is a Monte Carlo code for the modeling of protein sequence data 
(and more generally multivariate discrete data).
The basic usage of the code consists in 1) learning a probabilistic model 
for the joint distribution of variables from the correlations in real samples 
and 2) generating artificial discrete data sampled from the model via Markov chain Monte Carlo (MCMC) sampling.

**elss** was originally developed to show that pairwise models for protein sequences with 
correlated amino acids can be learned and resampled using MCMC methods 
([paper](http://www.pnas.org/content/112/44/13567)).

## Obtaining the source

All **elss** source code is hosted on Github. 
The latest version can be downloaded using 
[this link](https://github.com/simomarsili/elss/archive/v0.3.1.tar.gz). 

## Prerequisites

In order to compile **elss**, you will need a **Fortran compiler** installed on your machine.
If you are using Debian or a Debian derivative such as Ubuntu, you can install the gfortran
compiler using the following command:

```bash
$ sudo apt-get install gfortran
```

The inference algorithm works by simulating a swarm of persistent Markov chains. 
To compile **elss** with support for parallel runs on a distributed-memory architecture,
you will need to have a valid **MPI implementation** installed on your machine. 
The code has been tested and is known to work with the latest versions of both OpenMPI and MPICH.   

OpenMPI (recommended) can be installed on Debian derivatives with:
```bash
$ sudo apt-get install openmpi-bin libopenmpi-dev
```
For details on running MPI jobs with OpenMPI see [this link](https://www.open-mpi.org/faq/?category=running)

Alternatively, MPICH can be installed with:
```bash
$ sudo apt-get install mpich libmpich-dev
```

The compiling and linking of source files is handled by **Gnu Make**. 
If you are using Debian or a Debian derivative such as Ubuntu, you should
find Gnu Make 4.1 already installed. 

(optional) **git** version control software for obtaining the source code:
```bash
$ sudo apt-get install git
```

## Compiling and installing

To compile **elss**, type `make` in the `src` directory:
```bash
$ cd src; make
```
This will build the `elss` executables (`elss-learn`, `elss-sample` and `elss-eval`).

To install the executables, type `make install`:
```bash
$ cd src; make install
```
for the default installation dir (`/usr/local/bin`)
or use `DESTDIR` to override it.  
For example, to install in `~/.local` instead of `/usr/local`:
```bash
$ cd src; make install DESTDIR=~/.local
```

## Testing

Run the `run-test.bash` script in the test directory: 
```bash
$ cd test; bash run-test.bash
```

## Input data format

The input data should be encoded as space/tab separated integer labels,
with variables as columns and samples as rows:
```bash
$ head encoded.txt
20 13  8 21  4 17 18  8 19 24
12  0 13  0  6  4 12  4 13 19
19  4  2  7 13 14 11 14  6 24
 6 14 21  4 17 13 12  4 13 19
 3  4 15  0 17 19 12  4 13 19
 2  0 19  4  6 14 17  8  4 18
 2 14 13  3  8 19  8 14 13 18
```

The code assumes that all variables share a common set of classes, 
encoded as integer indices starting from zero.  
Alternatively, biological sequence data can be directly read from a
multiple sequence alignment file in FASTA format.

## Basic usage

The standard workflow has two steps:
1) the fitting of the pairwise model to data (using the `elss-learn` executable) and
2) the sampling of artificially generated data from the fitted model (using `elss-sample`).  

### elss-learn
The fitting consists in a first-order iterative minimization
of a cost function including two terms, a term proportional to the 
parameters' likelihood and a regularization term.
Example:
```bash
$ mpiexec -n 4 elss-learn --fasta 1.fa --niter 2000 -n 10000
```
- `mpiexec -n 4`: compute the gradient of the cost function simulating 4 independent Markov chains
- `--fasta 1.fa`: read data from file `1.fa` in FASTA format
- `--niter 2000`: set the number of iterations for iterative minimization to 2000
- `-n 10000`: set the length of each MC chain (per gradient evaluation) to 1000 MC sweeps

The run will produce a binary checkpoint file `chk`,
that contains all the fitted parameters.
For a full list of options and details, type `elss-learn -h`.

### elss-sample

`elss-sample` reads the parameters contained in the checkpoint file produced by `elss-learn` 
and simulate a MC trajectory sampling from a model of pairwise-interacting variables:
```bash
$ elss-sample --chk chk -n 100000 -u 100
```
- `--chk chk`: read the fitted model from checkpoint file `chk`
- `-n 100000`: run a MC trajectory for 100000 MC sweeps
- `-u 100`: dump a configuration every 100 MC sweeps to a `trj` file.  
The output of the calculation is a `trj` file containing 1000 (100000/100)
configurations sampled according to the fitted pairwise model.
For a full list of options and details, type `elss-sample -h`.

### The checkpoint (.chk) files
A checkpoint file is an unformatted binary file containing a set of fitted
parameters and all the informations needed to restart a
previously interrupted optimization. For example, the command:
```bash
$ mpiexec -n 4 elss-learn --chk old.chk --fasta 1.fa --niter 2000 -n 10000
```
will restart the optimization process from the values found in `old.chk`.

A checkpoint file can be converted to a plain text file using the `elss-pchk` tool together with the `-u` option.
The command `elss-pchk -u <file>` will generate a file named `<file>.txt`, that is:
```bash
$ elss-pchk -u chk
$ head chk.txt
        10 # n. of variables
        21 # n. of classes
   protein # n. output data format
         4 # n. samples in checkpoint
# samples start here
12  4  4  1 17  4  4  3 10  9
12 20  3 18 12  4  4  4 18  9
 7 14 16 18 17 13  9  3  5  3
 6 11  3 17 17  4  3 10 18  4
# parms start here
           1  0.36794536927277111      -0.38797251076999351       0.45155637053903591      ....
```
Viceversa, a custom text file containing user-defined parameters can be used 
to generate a checkpoint file using the `-f` option. Given a valid input file <file> of parameters, 
the command `elss-pchk -f <file>` will generate an unformatted checkpoint file named `<file>.chk`.

The lines of a valid input file to `elss-pchk -f` will contain, in this order:
- the number of features/variables in the system, NF
- the number of classes or possible values that can be taken by a variable, NC
- a keyword selecting the format of output data (int, protein)
- the number of samples contained in the subsequent lines, NS
- NS lines, each containing a sample encoded as a space/tab separated array of integer labels
- an arbitrary number of lines each containing the biases for each class of a given variable, with this format:  
  p x(1) x(2) ... x(NC)  
  where the {x} are the elements of the NC-long array of biases for variable p.  
  The program will set to zero the biases for those variables that are not explicitly defined.
- an arbitrary number of lines each containing the matrix of couplings for a pair of variables, with this format:  
  p q x(1,1) x(1,2) ... x(1,NC) ...  
  where the {x} are the elements of the NC x NC array of couplings for the variables p and q, iterated sequentally in row-major order.  
  **NB: q > p**  
  The program will set to zero the couplings for those pairs that are not explicitly defined.
- all characters following the `#` symbol are comments and are ignored
- empty lines are ignored

# Contributing

**elss** is an OPEN Source Project so please help out by [reporting bugs](http://github.com/simomarsili/elss/issues) or [forking and opening pull](https://github.com/simomarsili/elss) requests when possible.

# LICENSE (BSD 3 clause)

Copyright (c) 2016, Simone Marsili  
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Citation note

If you are going to publish results obtained using the software, please cite the **elss** papers:

[From residue coevolution to protein conformational ensembles and functional dynamics](http://www.pnas.org/content/112/44/13567)  
L. Sutto, S. Marsili, A. Valencia, F.L. Gervasio  
Proceedings of the National Academy of Sciences 112(44):13567-13572 (2015)
