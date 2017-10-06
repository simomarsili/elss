#!/bin/bash
tmp="$(mktemp)"
for ffile in *.f90;
do
    if head -1 "$ffile" | grep -q "Copyright (C)"; then
	echo $ffile;
	tail -n+4 $ffile > $tmp; cat HEADER.txt $tmp > $ffile
    fi
done
