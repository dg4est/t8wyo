!------------------------------------------------------------------------------
subroutine readpm(fileinput,              &
                  case_title,             &
                  mesh_file,              &
                  bcs_file,               &
                  fmach,yangle,zangle,re_number)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"

  character*200 fileinput
  character*200 case_title
  character*200 mesh_file
  character*200 bcs_file
  real(r8)    :: fmach,yangle,zangle,re_number
!--Tmps
  integer(i4) :: IFILE
  integer(i4) :: ierr
  integer(i4) :: ipass
  integer(i4) :: ilen
  logical     :: file_exists
      
!
!-------------------------SINGLE PE READS IN FILE-------------------------
!
      IFILE = 0
      if (ID_PROC .eq. 0) then
      write(iwrit,601) fileinput
  601 format(' Reading input list from file: ',a80)
!
      open(unit=iread,file=fileinput,form='formatted',status='old',err=800)
      IFILE = 1
      read(iread,525) case_title

      read(iread,500)
      read(iread,*) fmach
      read(iread,500)
      read(iread,*) yangle
      read(iread,500)
      read(iread,*) zangle
      read(iread,500)
      read(iread,*) re_number

      read(iread,500)
      read(iread,525) mesh_file
      call char_adjust(mesh_file)

      read(iread,500)
      read(iread,525) bcs_file
      call char_adjust(bcs_file)

  800 continue
      close(iread)

      endif
!--------------------------------END SINGLE PE READ-------------------------
!----------------------------------------------------------------------------
!
#ifdef MPI_ON
!
!--Now BroadCast All To other PEs
!
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
      call MPI_BCAST(IFILE    ,1     ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(case_title,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mesh_file, 200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bcs_file  ,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fmach     ,1  ,MPI_REAL8    ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(yangle    ,1  ,MPI_REAL8    ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(zangle    ,1  ,MPI_REAL8    ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(re_number ,1  ,MPI_REAL8    ,0,MPI_COMM_WORLD,ierr)

#endif
!----------------------------------------------------------------------
!
      IF (IFILE .ne. 1) THEN
        if (ID_PROC .eq. 0) THEN
          write(iwrit,611)
  611     format('INPUT FILE NOT FOUND ')
        ENDIF
        call stop_all
      ENDIF
!----------------------------------------------------------------------
      RETURN
!----------------------------------------------------------------------

  500 format(1x)
  525 format(a200)

end subroutine readpm
