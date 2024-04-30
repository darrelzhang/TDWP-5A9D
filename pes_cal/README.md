# PES calculation
This folder contains the program for calculating the potential energy surface of the X+YNZ reaction.

As shown in potdriv.f, the interface for the PES is given in the subroutine potdriv(...). 
The PES we used for H+NH3 is given in "nh4_pipnn.o". 
Due to copyright issues with the potential energy surface, we are unable to provide its source code.
Based on the requirements of the interface, one can change the PES for any reaction of the X+YNZ2 type.

# Requirements
Linux OS

Intel fortran compiler (2015 or later)

# compile
One should modify the variable MKL_PATH in the makefile and then

make

to obtain the executable file: xynh2_pes.exe

# run

Here we give a bash script "runpes.sh" for running.

In the scattering calculation, the potential energy surface is divided into the asymptotic region and the interaction region. These are represented by the "asy" and "int" options in the runpes.sh script, respectively. To calculate the potential energy surface faster, we use multiple CPUs to accelerate the computation, as indicated by the variable "potnum" in the runpes.sh script. One can change the value of this variable based on the number of available CPUs. Since the potential energy surface index starts from 0, “potnum” is equal to the nCPUs-1.

After the program runs, it will generate POT* potential files in the "result" folder. We need to copy these files to the "scatt/result" folder, as the scattering calculation will read these potential files.
