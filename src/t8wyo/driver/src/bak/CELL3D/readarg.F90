!----------------------------------------------------------------------
subroutine readarg(fileinput)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"

  character *200 fileinput

!Tmps...
  integer(i4) :: ISTOP,INPUT_VARS,ierr
  integer(i4) :: iargc

!-------------------------------------------------------------------------
#ifdef MPI_ON
!-------------------------------------------------------------------------
      ISTOP      = 0
      INPUT_VARS = 0   !Use for consistency with nsu3d for now...
!!    INPUT_VARS = 4   !For mpirun -np 4 ./pre_nsu3d  this should be the value (!!)
#if defined (mpxlf90_r) || (mpxlf90)
      INPUT_VARS = 0
#endif

      if (id_proc .eq. 0) then
      if ( iargc() .gt. INPUT_VARS ) then
        call getarg(1,fileinput)
      else
#if defined (mpxlf90_r) || (mpxlf90)
        write(iwrit,601)
#else
        write(iwrit,602)
#endif
  601   format(/' IBM(AIX) Parallel Mode Usage: poe nsu3d -procs ncpus [input_file]',/)
  602   format(/' Parallel Mode Usage: mpi -np ncpus nsu3d [input_file]',/)
        ISTOP      = 1
      endif
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ISTOP    ,1     ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fileinput,200 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      if (ISTOP .eq. 1) then
        call stop_all
      endif

      return
!-------------------------------------------------------------------------
#else
!-------------------------------------------------------------------------
      INPUT_VARS = 0
      if ( iargc() .gt. INPUT_VARS ) then
        call getarg(1,fileinput)
      else
        write(iwrit,601)
  601   format(/' Sequential Mode Usage: nsu3d [input_file]',/)
        stop
      endif

!-------------------------------------------------------------------------
#endif
!-------------------------------------------------------------------------
end subroutine readarg
