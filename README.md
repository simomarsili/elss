# **elss**

**elss** generates multivariate discrete data that preserve the correlations
observed in a real set of samples. The probabilistic model used for the generation
of artificial samples takes into account only pairwise interactions among variables.
Learning of the model's parameters is carried out using maximum a posteriori
(MAP) estimation and a stochastic MCMC approximation of the gradient of the cost
function.

<!---
**elss** has been recently used to show that a minimal model of
amino acids interacting in pairs is able to capture higher-order correlations
from sequences of proteins with a common ancestor (i.e. the presence of
clusters of sequences with specific biological functions, see our paper
[From residue coevolution to protein conformational ensembles and functional dynamics](http://www.pnas.org/content/112/44/13567)).
-->

## Obtaining the source

All **elss** source code is hosted on Github. 
You can download the latest version of the code using
[this link](https://github.com/simomarsili/elss/archive/v0.2.1.tar.gz). 

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

## Compiling

To compile **elss**, type `make` in the `src` directory:
```bash
$ cd src; make
```

This will build the `elss` executables (`elss-learn`, `elss-sample` and `elss-eval`).

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

The code assumes that all variables share a common set of classes.  
Aletrnatively, **elss** can read directly biological sequence data from a
multiple sequence alignment file in FASTA format:
```bash
$ head data.fa
>G3VGV0_SARHA/14-174
LAVMGTCCVGKTALTIQFTKNRFMTQYNPTCQDFYRKHTVADEERRQLDIVDTTSTEAFYCLRDQAMRWGEGFLLVYSVNDPHSFENVNVLWDHLQKLKGRVPMVLVANKVDVTDRLVNPRQGQEVARRFGVPYVETSAKSKEGVEQAFHELV
>A2FI73_TRIVA/11-193
VVTIGETAVGKTSIISRLVNARFSENESPTIGNFLMHEENIGNQKIELQIWDTAGQEKYRALSPIYCRDAAVGLIIYDVTNKDTFNKIDNWIKLFKDVADEALVYIVGNKCDKIELTVERNAIE-VFSDQGYNCFFTSAKTGEGINDLFHDIC
>T0LCH1_9MICR/12-173
IAILGYYSVGKSSLSLKYVRNQFNPNEESTIASYLTKSMSTKDSTIQFEIWDTAGQERYNSLVSIYYKNADAALIVYDITSRDSFEAAKQWVYELNFQKPDFLKILVGNKTDMEERQVDFEEGKEYAMQQNLIFLEASAKSGENVSKIFELFA
>Q3SDV0_PARTE/13-174
IVLVGDSGSGKTTLFMKHAEQQFCQNLSPTIIEFHNKFVEYQRKMIKLQLWDTAGQETFRSISQNYYRKANSIFFIYDITNKQSFERVYQWMNEAKQLAPDLIKVLIGNKSDLINRQVSFDEGKLFALENDLEFFELSAFGNRNLEDPIYYVL
>IFT27_TRYB2/7-170
VAVVGAPTVGKTAFVQMLHSNTFPKNYLMTLCDFIVKEVPVDDNTVEMIIFDVSGQREYEPMVSSYLQNTAVFIVMYDVSNKVTFEACARWVNQVRTNSK-SVGILIANKSDLSDAEVTDRQGKDLANANKMKFYKISTLRGVGITEPIDEIA
```

## Basic usage

The standard workflow starts by fitting a model to data (using `elss-learn`)
and then sampling artificial data from the fitted model (using the `elss-sample`
executable). We will discuss in detail the example calculations used as tests
(see `test` directory and `run-test.bash`). For a full list of options and 
more details, run the executables with the `-h` flag *e.g.* `elss-learn -h`.


### elss-learn
The fitting process corresponds to a first-order iterative minimization
of a cost function including two terms, one proportional to the
likelihood of the parameters and a regularization term.

```bash
$ mpiexec -n 4 elss-learn --fasta 10.fa --niter 2000 -n 1000
```
- `mpiexec -n 4`: compute the gradient using 4 independent Markov chains
- `--fasta PF00076.fa`: read data from file `PF00076.fa` in FASTA format
- `--niter 2000`: set the number of iterations
- `-n 1000`: set the length of each simulated trajectory to 1000 MC sweeps

The run will produce a binary checkpoint file `chk`,
that contains all the fitted parameters.

### elss-sample

The `chk` file can be used as input to `elss-sample`:
```bash
$ elss-sample --chk chk -n 100000
```
- `--chk chk` read the fitted model from checkpoint file `chk`
- `-n 100000`

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

