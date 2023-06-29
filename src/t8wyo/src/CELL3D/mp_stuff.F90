module mp_stuff

  use my_kinddefs

  real(r8)    ::  time_init,time_init_mp
  real(r8)    ::  time1,time1_mp
  real(r8)    ::  time2,time2_mp
  real(r8)    ::  time0(10),time0_mp(10)
  character*80::  ctime(10)
  integer(i4) ::  number_timings
  integer(i4) ::  MPI_BARRIER_LEVEL
  integer(i4) ::  id_proc,num_proc,io_debug1
  integer(i8) ::  ntotal_dyn_mem

end module mp_stuff

