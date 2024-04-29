objdir:=obj
moddir:=mod
srcdir:=src

FC      :=mpif90
modfile :=-module ./$(moddir)
FFLAGS  :=-O2 -132 -traceback -openmp  $(modfile) 
MKL_PATH:=/work1/soft/intel2015/composerxe/mkl/lib/intel64
LDFLAGS :=-L$(MKL_PATH)  -lmkl_rt


pes0:=nh4_pipnn.o 


objs=$(addprefix $(objdir)/, varmd.o readinp.o gaussq.o lib.o rotbas.o abs_flux.o \
	         array_zero.o trans_sin_dvr.o small.o  ${pes0} \
                 ctrans.o degree_order.o exp_eig.o init_wp.o potdriv.o \
                 nh2_hamiltonian.o nh3_ref_pot.o nh3_ref_ham.o nh3_hamiltonian.o \
                 nh3_grid_matrix.o nh4_rotmatrix.o rotmat_vpot.o trans1.o \
                 scatt_nh2.o plotwave.o schprog.o reaction.o strans0.o  \
                 mpi_potent.o mpi_scatt.o main.o )

xynh2.exe: $(objs)
	 $(FC) $(FFLAGS) $^ $(LDFLAGS) -o $@

$(objdir)/%.o: $(srcdir)/%.f
	$(FC) $(FFLAGS) -c $< -o $@

$(objdir)/%.o: $(srcdir)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@


clean:
	@mv obj/nh4_pipnn.o obj/nh4_pipnn.o.bak
	rm -f $(moddir)/*.mod $(objdir)/*.o *.exe
	@mv obj/nh4_pipnn.o.bak obj/nh4_pipnn.o





