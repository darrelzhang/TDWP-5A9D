objdir:=obj
moddir:=mod
srcdir:=src

exefile:=ynh2bound.exe

FC      :=ifort
pes_mod :=-module ./mod
FFLAGS  :=-O2 -132 -traceback  -openmp  $(pes_mod) 
MKL_PATH:=/work1/soft/intel2015/composerxe/mkl/lib/intel64
LDFLAGS := -L$(MKL_PATH)  -lmkl_rt 
LIBA    :=libarpack_intel64.a 

pes0:=nh4_pipnn.o

objs:=$(addprefix $(objdir)/, varmd.o small.o interpolation.o gaussq.o trans_sin_dvr.o \
                        integral.o $(pes0) potdriv.o nh2_hamiltonian.o \
                        nh3_ref_pot.o nh3_ref_ham.o nh3_hamiltonian.o \
                        nh3_grid_matrix.o init.o array_zero.o \
                        trans0.o trans1.o nh3_rotmatrix.o \
                        hmatdriv.o arpack.o main.o )

$(exefile): $(objs)
	$(FC) $(FFLAGS) $^ $(LDFLAGS) $(LIBA) -o $@

$(objdir)/%.o:  $(srcdir)/%.f
	$(FC) $(FFLAGS) -c  $< -o $@
$(objdir)/%.o:  $(srcdir)/%.f90
	$(FC) $(FFLAGS) -c  $< -o $@

clean:
	@mv $(objdir)/nh4_pipnn.o $(objdir)/nh4_pipnn.bak
	rm -f  $(moddir)/*.mod $(objdir)/*.o $(exefile)
	@mv $(objdir)/nh4_pipnn.bak $(objdir)/nh4_pipnn.o 
