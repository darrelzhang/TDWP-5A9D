
# PES calculation
This folder contains the program for calculating the potential energy surface of the X+YNZ reaction.

As shown in potdriv.f, the interface for the PES is given in the subroutine potdriv(...). 
The PES we used for H+NH3 is given in "nh4_pipnn.o". 
Due to copyright issues with the potential energy surface, we are unable to provide its source code.
Based on the requirements of the interface, one can change the PES for other type of X+YNZ2 reaction.

# Requirements
Linux OS

Intel fortran compiler (2015 or later)

# compile
One should modify the variable MKL_PATH in the makefile and then

make

to obtain the executable file: xynh2_pes.exe

# run

mpirun -hostfile machinefile -n  procs  ./xynh2.exe

machinefile： contains the hostname of the parallel computing nodes
