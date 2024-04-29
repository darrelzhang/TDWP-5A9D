objdir:=obj
moddir:=mod
srcdir:=src

FC      :=ifort 
modfile :=-module ./$(moddir)
FFLAGS  :=-O2 -132 -traceback -openmp  $(modfile) 
MKL_PATH:=/work1/soft/intel2015/composerxe/mkl/lib/intel64
LDFLAGS :=-L$(MKL_PATH)  -lmkl_rt

pes0:=nh4_pipnn.o

objs=$(addprefix $(objdir)/, varmd.o readinp.o gaussq.o lib.o rotbas.o \
	         array_zero.o trans_sin_dvr.o small.o  $(pes0) \
                 degree_order.o init_wp.o potdriv.o \
                 nh2_hamiltonian.o nh3_ref_pot.o nh3_ref_ham.o main.o )
                   

xynh2_pes.exe: $(objs)
	 $(FC) $(FFLAGS) $^ $(LDFLAGS) -o $@

$(objdir)/%.o: $(srcdir)/%.f
	 $(FC) $(FFLAGS) -c $< -o $@

$(objdir)/%.o: $(srcdir)/%.f90
	 $(FC) $(FFLAGS) -c $< -o $@

clean:
	@mv obj/nh4_pipnn.o obj/nh4_pipnn.o.bak
	rm -f $(moddir)/*.mod $(objdir)/*.o *.exe
	@mv obj/nh4_pipnn.o.bak obj/nh4_pipnn.o
