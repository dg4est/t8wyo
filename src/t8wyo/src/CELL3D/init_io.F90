!-----------------------------------------------------------------------------
subroutine init_io

   
  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none

  integer(i4) :: k

!--Version

  IO_VERSION = 3

!--Level of reporting  (0,1,2)

  IVERBOSE   = 1
!
!--Level of MPI-BARRIER (high = more barriers)

  MPI_BARRIER_LEVEL = 3

!--Parallel I/O Settings
!
  IPARALLEL_READ  = 1
  IPARALLEL_WRITE = 1
  IPARALLEL_AMG_WRITE = 1
!
!--Set IO_SYSTEMCALL to control std output messages
!
  IO_SYSTEMCALL = 2 !Possible Values: 0,1,2

!
!--Wait Time for NFS to update Cache after id_proc=0 makes directory for IO (secs)
!
  IO_SLEEP      = 0

!
!--Set I/O Numbers
!
  do k=1,100
   IO_UNITS(k) = k + 100
  enddo

! io_mem = 6
  io_mem = 128
  io_mem_proc = 1

  io_mesh  = 101
  io_mcell = 21
  io_poin1 = 22
  io_distf = 23
  io_bcs   = 24
  io_comp  = 25
  io_amg   = 26
  io_local = 200

end subroutine init_io
