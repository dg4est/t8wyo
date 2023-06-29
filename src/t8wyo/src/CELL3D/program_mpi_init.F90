!------------------------------------------------------------------------------
subroutine program_mpi_init
!
!--Initialize MPI and I/O Unit numbers
!--report on mpi initialization
!--and initialize timing   
!--NOTE:  THIS MUST BE THE FIRST EXECUTABLE CALLED BY MAIN PROGRAM....
!
  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"

 
  integer(i4) :: ierr
 
#ifdef MPI_ON
      !call MPI_Init(ierr) ! AK: initialized in t8wyo
      call MPI_Comm_rank(MPI_COMM_WORLD, id_proc , ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierr)
#else
      id_proc  = 0
      num_proc = 1
#endif

     ntotal_dyn_mem = 0
 
!------------------------------------
!--Set BASIC I/O Unit Numbers
!
!CDX1*iread    = 5
      iread    = 55
      iunbuf   = 0
      iwrit    = 6
      iterm1   = 5
      iterm2   = 6
!------------------------------------

    if (id_proc .eq. 0) then
    write(iwrit,601)
601 format(//'------------------------------------------------------------------------------',&
           //'                 ',                                                             &
             ' *** PARALLEL PROCESSING SETUP INFORMATION *** ',                               &
            /'                     ',/)
    end if

#ifdef MPI_ON
    if (id_proc .eq. 0) then
    write(iwrit,602) num_proc
602 format(//'  MPI RUN INITIATED; No of Processors: ',i7,//)
    end if
#endif

    number_timings = 0
    call get_time_mp(time_init,time_init_mp)

end subroutine program_mpi_init
!------------------------------------------------------------------------------
