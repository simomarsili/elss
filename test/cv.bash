#!/usr/bin/env bash 
# compute CV vs T given a rst file

n=10
EXE=../src/elss
RST=rst
for i in {1..40..1};
do
    t=$(echo "$i/20.0"|bc -l);
    $EXE -r rst --nsweeps 1000000 -t $t;
    # print t, <E>/N, Cv/(N*Kb)
    grep ">" 0.trj | awk -v a=$t -v n=39 '{e=$3;m+=e;m2+=e**2}END{m=m/NR;m2=m2/NR;print a, m/n, (m2-m**2)/(n*a**2)}' ;
done 

