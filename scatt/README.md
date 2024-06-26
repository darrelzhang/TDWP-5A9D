# TDWP-5A9D
a parallel algorithm for high-dimensional quantum dynamics  simulations in poly-atomic reactions

The parallelization strategy for wavepacket dynamics calculations studying polyatomic reactions is given in the program mpi_scatt_1.f and mpi_scatt_2.f.

As shown in potdriv.f, the interface for the PES is given in the subroutine potdriv(...). 
The PES we used for H+NH3 is given in "nh4_pipnn.o". 
Due to copyright issues with the potential energy surface, we are unable to provide its source code.
Based on the requirements of the interface, one can change the PES for any reaction of X+YNZ2 type.

# Requirements
Linux OS

Intel fortran compiler (2015 or later)

openmpi-1.4.5 (or later)

*We carried out calculations on one of our in-house computer clusters. The cluster comprises 30 shared-memory computer nodes, each with 20 Intel Xeon E5-2640 v4 processors (CPU cores) at 2.40 GHz and 256 gigabytes (GB) of main memory. The operating system running on the nodes is CentOS Linux release 7.2.1511 with the 3.10.0-327.el7.x86_64 kernel. Each node is equipped with one Infiniband Host Channel Adapter (HCA) supporting 4x Quad Data Rate (QDR) connections with a 40 GB/second speed. We have one 4x-QDR Infiniband connection from each node to a central Infiniband switch (Mellanox, MT26428). The TDWP-5A9D program was compiled at the -O2 optimization level with OpenMP support using the Intel Fortran compiler (version 15.0.2). The MPI utilities implemented as in Open MPI 34 (version 1.4.5) were adopted.*


# compile
One should modify the variable MKL_PATH in the makefile and then

make

to obtain the executable file: xynh2.exe

# run

mpirun -hostfile machinefile -n  procs  ./xynh2.exe

machinefile： contains the hostname of the parallel computing nodes

procs: number of MPI processes

We also give a bash script "run.sh" in the folder. The run options "asy", "int", and "sca" represent reading the asymptotic region, the interaction region potential energy surface, and the time-dependent wave packet scattering，respectively.
In "input" file, the parameter "npot" is equal to the number of "POT*" files (i.e.ncpus. See the README in "pes_cal" folder.)

The input file for the program is given as "input" in default.
All of the output files are in "result" directory.
