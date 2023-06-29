!-------------------------------------------------------------------------------
subroutine abort_all

  use my_kinddefs
  use mp_stuff
  implicit none

#include "mympif.h"

  integer(i4) :: ierrocode,ierr

          ierr = 0
          ierr = 1/ierr
#ifdef MPI_ON
   call MPI_Abort(MPI_COMM_WORLD, ierrocode,ierr)
   call MPI_Finalize(ierr)
#endif

   stop

end subroutine abort_all
