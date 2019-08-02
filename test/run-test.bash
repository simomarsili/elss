#!/usr/bin/env bash

HEADER=$\
'
======================
Requirements/Tested on
======================

- GNU Make 3.81 or higher (4.1)
- gfortran 4.8.2 or higher (7.4.0)
- OPENMPI 1.6 or higher (2.1.1)
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

# always compile
(cd $src_dir &&
	command -v mpif90 >/dev/null 2>&1 || { echo >&2 "mpif90 is required but it's not installed.  Aborting."; exit 1; }
 echo "compiling elss..."
 make realclean; 
 make; 
) &> log ;

EXE=$src_dir/elss

echo '
================
Running tests...
================
'

rm chk trj ene sample.test.out eval.test.out

command -v mpiexec >/dev/null 2>&1 || { echo >&2 "mpiexec is required but it's not installed.  Aborting."; exit 1; }

echo "fitting model parameters... (chk, prm, log)"
mpiexec -n $NPROC $EXE\-learn --fasta $DATA -n $NS --seed 123 >> log 2>&1;
mv log learn.log;

echo "sampling sequences from the model distribution... (trj, log)"
$EXE\-sample -c chk -n 100000 --seed 123 >> log 2>&1;
tail -20 trj > sample.test.out;
mv log sample.log;

echo "checking data energies... (ene, log)"
$EXE\-eval -c chk --fasta $DATA >> log 2>&1;
tail -10 ene > eval.test.out;
mv log eval.log;

#echo "checking prms... "
#$EXE\-pchk -u chk > chk.dat; 

if [ $NITER -eq 2000 ] && [ $NS -eq 1000 ] && [ $NPROC -eq 4 ] && [ $DATA == '10.fa' ]
then 
    echo '
================
Checking results
================
'
    # check diffs
    tests_ok=true
    
    if ! cmp sample.test.out sample.ref.out >/dev/null 2>&1; then
	echo "elss-sample: RESULTS DIFFER..."
	echo "compare sample.test.out and sample.ref.out for minor numerical diffs and sample.log file"
	tests_ok=false
    fi
    
    if ! cmp eval.test.out eval.ref.out >/dev/null 2>&1; then
	echo "elss-eval: RESULTS DIFFER..."
	echo "compare eval.test.out and eval.ref.out for minor numerical diffs and eval.log file"
	tests_ok=false
    fi

#    if ! cmp chk.dat CHK.dat >/dev/null 2>&1; then
#	echo "ene: RESULTS DIFFER..."
#	echo "check files chk.dat and CHK.dat for minor numerical diffs and log file"
#	tests_ok=false
#    fi

    if [ "$tests_ok" = true ] ; then
	echo "!!! TESTS OK !!!"
	rm sample.test.out eval.test.out
    fi
fi 

exit 0
