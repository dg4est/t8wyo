!-------------------------------------------------------------------------------
subroutine stop_all

  use my_kinddefs
  implicit none
#include "mympif.h"

  integer(i4) :: ierr

#ifdef MPI_ON
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Finalize(ierr)
#endif

  stop
  end
