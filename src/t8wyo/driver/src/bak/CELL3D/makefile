#------------------------------------------------------------------------------
#--->SET UTILITY NAMES (machine dependent)

#NOTE: -O0 compiler optimiation required for: build_tran*.f90 and modify_real_ghost_transfer_mpi.f90
#for older versions of ifort

#NOTE: -O0 compiler optimization required for: build_tran*.f90 and modify_real_ghost_transfer_mpi.f90
#for older versions of ifort

##MPI COMPILERS (Set appropriate flags for underlying compiler below)
FC     = mpiifort
CC     = mpiicc
#FC     = mpif90
#CC     = mpicc


# Intel Compiler Flags
#FFLAGS =  -r8 -O0 -convert big_endian
FFLAGS =  -r8 -O0 -g -check all -fpe0 -traceback -convert big_endian  -assume buffered_io
FFLAGS_NOOPT = -r8 -O0 -fpe0 -traceback -convert big_endian
#FFLAGS =  -check all -r8  -O0 -mieee-fp -fpe0 -traceback -w  -convert big_endian

# GCC Fortran Compiler Flags
#FFLAGS = -O3 -march=native -fdefault-real-8 -fdefault-double-8 -fconvert=big-endian -ffree-line-length-none #-frounding-math -fsignaling-nans

# CPP Options
CPP    =  cpp -traditional -ffreestanding
#CPP = cpp -C -P -traditional

#Note that INCLUDEPATH is always set to Metis directory since this specifies
#the Metis include file location for the C-Wrapper code only
INCLUDEPATH = -I$(PARMETISDIR)/include
#Set options for including partitioning libraries
#INCLUDELIB  = -L$(PARMETISDIR)/lib -lparmetis -lmetis -lGKlib -L../drums/lib -ldrums  #METIS GIT Version (for gfortran)
#INCLUDELIB  = -L$(PARMETISDIR)/lib -lparmetis -lmetis  -L../drums/lib -ldrums #V4.0.3

 

#------------------------------------------------------------------------------
#SCRIPT CONFIGURABLE COMPILER
#FC#
#FFLAGS#
#CPP#
#------------------------------------------------------------------------------
#--->SET OPTIONS

#ARCH : Choose from: SV1, T3E, IBM, INTEL, SGI, OTHER
ARCH = INTEL

#VEC  : Choose from: VEC_ON, VEC_OFF
VEC  = VEC_OFF

#MPI  : Choose from: MPI_ON, MPI_OFF
MPI  = MPI_ON

#MPITHREAD : Choose from MPITHREAD_ON MPITHREAD_OFF
MPITHREAD = MPITHREAD_OFF

#SAFE : Choose from SAFE_ON SAFE_OFF
SAFE=SAFE_ON



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#PROGRAM NAME: (for tarfile)
PNAME = cell3d.v1.0
#------------------------------------------------------------------------------

SRC= .src
CNTRL= .cntrl
SHELL =/bin/sh


OBJECT_MOD=  my_kinddefs.o mpi_schedule_typedef.o my_typedefs.o mpi_schedules.o io_params.o  mp_stuff.o \
             my_mpi_barrier.o params.o  get_filename1.o global_arrays_dynamic.o local_arrays_dynamic.o my_allocate.o mpi_schedule_read.o my_deallocate.o all_reduce_max.o all_reduce_min.o all_reduce_sum.o mpi_to_gpt.o triarea.o

OBJECT_F90= cell3d.o cell3d_initialize.o cell3d_init.o program_mpi_init.o stop_all.o abort_all.o get_time_mp.o readarg.o init_io.o readpm.o char_adjust.o wout.o cell3d_read_mesh.o cell3d_flux_test.o cell3d_test_connectivity.o 


cell3d  : $(OBJECT_MOD) $(OBJECT_F90)
	$(FC) $(FFLAGS) -o cell3d $(OBJECT_MOD) $(OBJECT_F90) $(INCLUDELIB)
	touch $(SRC)/empty


%.o: %.f90
	sed 's/#MPI/#/g' $*.f90 > $*.ftmp1
	sed 's/#IBL/#/g' $*.ftmp1 > $*.ftmp2
	$(CPP) -D$(FC) -D$(ARCH) -D$(MPI) -D$(MPITHREAD) -D$(VEC)  $*.ftmp2 > $(SRC)/$*.f90
	$(FC) -c $(FFLAGS)  $(ACT_DISK_MODULE) $(SRC)/$*.f90
	rm $*.ftmp1
	rm $*.ftmp2

%.o : %.c
	$(CC) -c $(CFLAGS) $(INCLUDEPATH) $(INCLUDELIB) $*.c


tarfile:
	@if test -d $(PNAME); then \
	continue; \
	else \
	mkdir $(PNAME); \
	fi
	cp -a *.f90 *.c *.h *.comm makefile README README.input input.pre_nsu3d.*  $(PNAME)
	mkdir $(PNAME)/.src
	tar cfvz $(PNAME).tgz $(PNAME)
	/bin/rm -rf $(PNAME)
	@echo "----------------------------------------"
	@echo "source tarfile is called $(PNAME).tgz"
	@echo "----------------------------------------"

clean:
	rm -rf $(SRC)/* *.o *.mod cell3d

