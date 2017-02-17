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
LAMBDA=0.01
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
	echo "compiling mcsg..."
	make realclean; 
	make; 
    fi
) &> log ;

EXE=$src_dir/mcsg

echo '
================
Running tests...
================
'

rm rst 0.trj 0.ene

command -v mpiexec >/dev/null 2>&1 || { echo >&2 "mpiexec is required but it's not installed.  Aborting."; exit 1; }

echo "estimating model parameters... (dump: rst,prm,LEARN.log)"
mpiexec -n $NPROC $EXE\-learn --fasta $DATA -n $NS --random_seed 123 --prefix learn >> log 2>&1;

echo "sampling sequences from the model distribution...(dumps: 0.trj,SIM.log)"
$EXE\-sample -r rst -n 100000 --random_seed 123 --prefix sample >> log 2>&1; 

echo "checking data energies...(dumps: 0.ene,EVAL.log)"
$EXE\-eval -r rst --fasta $DATA --prefix eval >> log 2>&1; 

if [ $NITER -eq 2000 ] && [ $NS -eq 1000 ] && [ $NPROC -eq 4 ] && [ $DATA == '10.fa' ]
then 
    echo '
================
Checking results
================
'
    # check diffs
    tests_ok=true
    
    if ! cmp 0.trj 0.TRJ >/dev/null 2>&1; then
	echo "0.trj: RESULTS DIFFER..."
	echo "check files 0.TRJ and 0.trj for minor numerical diffs and log file"
	tests_ok=false
    fi
    
    if ! cmp 0.ene 0.ENE >/dev/null 2>&1; then
	echo "0.ene: RESULTS DIFFER..."
	echo "check files 0.ENE and 0.ene for minor numerical diffs and log file"
	tests_ok=false
    fi

    if [ "$tests_ok" = true ] ; then
	echo "!!! TESTS OK !!!"
    fi
fi 

exit 0
