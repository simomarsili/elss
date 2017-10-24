# **elss**

**elss** can generate artificial discrete data preserving the correlations
observed in real data. The probabilistic model used for the generation of
artificial samples takes into account only pairwise interactions among variables.
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
The latest version can be downloaded using 
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

The code assumes that all variables share a common set of classes.  
Alternatively, biological sequence data can be directly read from a
multiple sequence alignment file in FASTA format.

## Basic usage

The standard workflow has two steps:
1) the fitting of the pairwise model to data (using the `elss-learn` executable) and
2) the sampling of artificially generated data from the fitted model (using `elss-sample`).  

### elss-learn
The fitting consists in the first-order iterative minimization
of a cost function including two terms, one proportional to the
likelihood of the parameters and a regularization term.
The `elss-learn` options control the parameters of the minimization run.
Example:
```bash
$ mpiexec -n 4 elss-learn --fasta 1.fa --niter 2000 -n 10000
```
- `mpiexec -n 4`: compute the gradient using 4 independent Markov chains
- `--fasta 1.fa`: read data from file `1.fa` in FASTA format
- `--niter 2000`: set the number of iterations to 2000
- `-n 10000`: set the length of each MC chain to 1000 MC sweeps

The run will produce a binary checkpoint file `chk`,
that contains all the fitted parameters.
For a full list of options and details, type `elss-learn -h`.

### elss-sample

The `chk` file can be used as input to `elss-sample`:
```bash
$ elss-sample --chk chk -n 100000 -u 100
```
- `--chk chk`: read the fitted model from checkpoint file `chk`
- `-n 100000`: run a MC trajectory for 100000 MC sweeps
- `-u 100`: dump a configuration every 100 MC sweeps to a `trj` file.
The output of the calculation is a `trj` file containing 1000 (100000/100)
configurations sampled according to the fitted pairwise model.
For a full list of options and details, type `elss-sample -h`.

### The checkpoint (chk) files
A `chk` file is an unformatted binary file containing the values for the
parameters of a fitted model and all the informations needed to restart a
previously interrupted optimization. For example, the command:
```bash
$ mpiexec -n 4 elss-learn --chk old.chk --fasta 1.fa --niter 2000 -n 10000
```
will restart the optimization process from the parameters' values found in `old.chk`.

A checkpoint file can be printed in a plain formatted ascii file using the `elss-pchk` tool together with the `-u` option:
```bash
$ elss-pchk -u chk > chk.txt
$ head chk.txt
# <data_type> <nvars> <nclasses> <ndata>
protein   10   21    4
 12  4  4  1 17  4  4  3 10  9
 12 20  3 18 12  4  4  4 18  9
  7 14 16 18 17 13  9  3  5  3
  6 11  3 17 17  4  3 10 18  4
           1  0.36794536927277111      -0.38797251076999351       0.45155637053903591        1.7133437108309755E-002 -0.41269123360182003       0.26054877643905250        7.6410686025211574E-002 -0.46495283238231327       -7.5388894712554749E-002 -0.26968317076052950      -0.34067001518143014       0.28254938511689753        1.2635022603526687      -0.21629051023435544      -0.22067848866868223       0.54605377753264595       0.22482120649044776      -0.33967234525728579      -0.50416522080998660      -0.34548697730341321        9.8838154756288280E-002
           2   3.8201087804195642E-002 -0.52011686577981775        2.5477078250094931E-002  0.31038273866588717       0.28780777424884108      -0.28573273645368163       -7.8576279813225833E-002 -0.34712250130692951       0.33571771352534041       -4.9860515476613354E-002 -0.45099425923628317       -3.6957228337975045E-002  0.55828017155671006        2.3965244789650966E-002   8.2699732321155364E-002  0.16096665766007959       -5.3525710910287357E-002 -0.24437923526905250       -2.1766753353324640E-002  0.30570751917965777       -7.4808760619281281E-002
           3   4.6460897226781803E-002 -0.45355018583070683       0.64617967902436269       0.51833678813212147      -0.28273631742632710       0.28794233802974700       -6.2872366540245056E-002 -0.38250247932398707       0.19011988646077030      -0.30648285926600810      -0.36951166191576035       0.41879206616547260      -0.45005440651910955        7.9741656538793407E-002  0.17473790680835091       0.44475601901850909       0.15913110729620758      -0.21501043124135605      -0.26056449326685033       -9.2738948099961041E-002 -0.13525447115267641     
           4  0.50677836706167712       0.12249587883667454      -0.42238086817404413      -0.45688575015404492        3.1246366771923467E-002 -0.26858330491332350      -0.27733409403105297       0.69788853345482527      -0.20276733289634646       0.48003547762691817       0.28028920246510725      -0.48738660590857491      -0.40806984768095667      -0.26550772980377579      -0.43633936054348721      -0.14070731240766471       0.76289426572635233       0.87537935147411239      -0.16659028270149304      -0.10396616229137706       -8.9761377651374735E-002
```
that will print to standard output a human-readable version of the actual content of file `chk`.  
A formatted file of parameters can be modified and a new unformatted checkpoint file can be generated using the `-f` option:
```bash
$ elss-pchk -f chk.dat > new.chk
```
A valid input to `elss-pchk -f` must contain:
- a line with
  - the data type (`protein`, `nuc_acid`, `int`)
  - the number of variables of the system
  - the number of possible values that a variable can take
  - the number of data samples (`ndata`) in the checkpoint
  in this order.
- a number `ndata` of lines containing different samples encoded as integer classes,  
  that will be used as starting points for the MC trajectories
- an arbitrary number of lines 
  






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

