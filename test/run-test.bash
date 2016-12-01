#!/usr/bin/env bash

HEADER=$\
'
======================
Requirements/Tested on
======================

- GNU Make 3.81 or higher (4.1)
- gfortran 4.8.2 or higher (5.3.1)
- OPENMPI 1.6 or higher (1.10)
'

DATA=10.fa
NS=1000
NITER=2000
NPROC=4

# check 
echo "$HEADER"

root_dir=../
src_dir=$root_dir/src

# check src dir
if [ ! -d $src_dir ]
then
    echo "run this script in test dir. Aborting."; exit 1; 
fi

# check executable 
(cd $src_dir && 
    if [ ! -f elss ]; then
	command -v mpif90 >/dev/null 2>&1 || { echo >&2 "mpif90 is required but it's not installed.  Aborting."; exit 1; }
	echo "compiling mcDCA..."
	make realclean; 
	make; 
    fi
) &> log ; 

EXE=$src_dir/mcdca

echo '
================
Running tests...
================
'

command -v mpiexec >/dev/null 2>&1 || { echo >&2 "mpiexec is required but it's not installed.  Aborting."; exit 1; }

echo "estimating model parameters... (dump: rst,prm,LEARN.log)"
mpiexec -n $NPROC $EXE --nsweeps $NS --fasta $DATA --learn-agd $NITER --random_seed 123 --lambda 0.01>> log 2>&1; 

echo "sampling sequences from the model distribution...(dumps: 0.trj,SIM.log)"
$EXE -r rst --nsweeps 100000 --random_seed 123 >> log 2>&1; 

echo "checking data energies...(dumps: 0.ene,EVAL.log)"
$EXE\-eval -r rst --fasta $DATA --random_seed 123 >> log 2>&1; 

if [ $NITER -eq 2000 ] && [ $NS -eq 1000 ] && [ $NPROC -eq 4 ] && [ $DATA == '10.fa' ]
then 
echo '
================
Checking results
================
'
    # check diffs 
    if ! cmp 0.trj 0.TRJ >/dev/null 2>&1; then
	echo "RESULTS DIFFER..."
	echo "check files 0.TRJ and 0.trj for minor numerical diffs"
	echo "check log file"
    else
	echo "!!! TESTS OK !!!"
    fi
fi 

exit 0
