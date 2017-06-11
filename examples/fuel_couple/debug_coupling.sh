#!/bin/bash

mpirun -np 1 xterm -e gdb --args  ../../fermi fuel_1.fer -c couple.dat : -np 1 xterm -e gdb control
