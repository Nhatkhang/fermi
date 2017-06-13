#!/bin/bash

# STEADY
#break_fer=( 'fer_couple.c:50' 'fer_main.c:25' ) 

# TRANSIENT
#break_fer=( 'fer_couple.c:181' 'fer_couple.c:50' 'fer_main.c:65' ) 
#break_fer=( 'fer_main.c:65' ) 
break_fer=( ) 

#break_con=( 'control.c:87' 'control.c:122' ) 
break_con=( 'control.c:123' ) 


# BREAKPOINTS FERMI
for i in ${break_fer[@]}
do
  gdbopt_fer="$gdbopt_fer -ex 'break $i' "
done
# auto run
gdbopt_fer+="-ex 'r' "

# BREAKPOINTS CONTROL
for i in ${break_con[@]}
do
  gdbopt_con="$gdbopt_con -ex 'break $i' "
done
# auto run
gdbopt_con+="-ex 'r' "

#debug command
if [ $# -eq 0 ];then
eval mpirun -np 1 xterm -e gdb $gdbopt_fer --args  ../../fermi iaea_coup.fer -c couple.dat : \
            -np 1 xterm -e gdb $gdbopt_con control
fi

#execution command
if [ ${#} -gt 0 ];then
mpirun -np 1 ../../fermi iaea_coup.fer -c couple.dat : \
       -np 1 control
fi
