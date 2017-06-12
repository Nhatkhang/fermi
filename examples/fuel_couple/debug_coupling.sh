#!/bin/bash

bps=( 'fer_couple.c:50' 'fer_couple.c:51') 


# BREAKPOINTS
for i in ${bps[@]}
do
  exopt="$exopt -ex 'break $i' "
done

echo $exopt
gdbcomm="gdb $exopt --args  ../../fermi fuel_1.fer -c couple.dat"
echo $gdbcomm

mpirun -np 1 xterm -e "$gdbcomm" : -np 1 xterm -e gdb control
