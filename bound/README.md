
# bound state for YNZ2
This folder contains the program for solving the bound states of YNH2 molecules.

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

to obtain the executable file: ynh2bound.exe

# run

./ynh2bound.exe
One can adjust the parameters in the "input" file for the six degree of freedoms.

The input file for the program is given as "input" in default.
All of the output files are in "result" directory.

For the "input" file, "nr" is the number of sine basis functions. "nrb" and "nvb" are the number of PODVR functions.

After the program runs, it will generate a "bound_wave" file, which stores the bound state wave function of the YNZ2 molecule. The corresponding energy levels are stored in the "bound_eigval" file. Subsequently, we need to copy the "bound_wave" file to the "scatt" folder, as it needs to be accessed during the scattering calculation.
