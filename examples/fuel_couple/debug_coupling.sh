#!/bin/bash

#set -x
# BREAKPOINTS
bps=( "fer_couple.c:50" ) 


# BREAKPOINTS COMMANDS
exopt=
for i in "${bps[@]}"
do
exopt="$exopt -ex 'break $i'"
done

# PRETTY OPTION
#exopt="$exopt -ex 'set print pretty on'"

echo $exopt
gdbcomm="gdb "$exopt" --args  ../../fermi fuel_1.fer -c couple.dat"
echo $gdbcomm

#set +x

#mpirun -np 1 xterm -e $gdbcomm : -np 1 xterm -e gdb control

mpirun -np 1 xterm -e gdb -ex 'break fer_couple.c:50'                    \
                          -ex 'set print pretty on'                      \
			  --args  ../../fermi fuel_1.fer -c couple.dat : \
       -np 1 xterm -e gdb -ex 'break control.c:105' control

#a="-ex 'break fer_couple.c:50'"
#mpirun -np 1 xterm -e gdb $a                                             \
#                          -ex 'set print pretty on'                      \
#			  --args  ../../fermi fuel_1.fer -c couple.dat : \
#       -np 1 xterm -e gdb -ex 'break control.c:105' control
