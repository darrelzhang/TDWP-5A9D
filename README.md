# TDWP-5A9D
a parallel algorithm for high-dimensional quantum dynamics  simulations in poly-atomic reactions

The parallelization strategy for wavepacket dynamics calculations studying polyatomic reactions is given in the program mpi_potent.f and mpi_scatt.f.

# Requirements
Linux OS

Intel fortran compiler (2015 or later)

openmpi-1.4.5 (or later)


# compile
One should modify the variable MKL_PATH in the makefile and then

make

# run

mpirun -hostfile machinefile -n  procs  ./xynh2.exe

machinefile： contains the hostname of the parallel computing nodes

procs: number of MPI processes

The input file for the program is given as "input" in default.
All of the output files are in "result" directory.
