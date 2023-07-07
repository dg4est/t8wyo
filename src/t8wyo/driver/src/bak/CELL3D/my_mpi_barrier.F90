!-------------------------------------------------------------------------------
module my_mpi_barrier_mod
contains

subroutine my_mpi_barrier(ilevel)

  use my_kinddefs
  use mp_stuff
  implicit none

  integer(i4), intent(in) ::  ilevel

#include "mympif.h"

  integer(i4) :: ierr

#ifdef MPI_ON
    if (ilevel <= MPI_BARRIER_LEVEL) then
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    endif
#endif


end subroutine my_mpi_barrier
end module my_mpi_barrier_mod
