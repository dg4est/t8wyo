module all_reduce_sum_mod

  interface all_reduce_sum

     subroutine all_reduce_sum_i4(iarray1,iarray2,ldim)
        use my_kinddefs
        implicit none
        integer(i4), intent(in)  :: iarray1(ldim)
        integer(i4), intent(out) :: iarray2(ldim)
        integer(i4), intent(in)  :: ldim
     end subroutine all_reduce_sum_i4

     subroutine all_reduce_sum_i41(iarray1,iarray2)
        use my_kinddefs
        implicit none
        integer(i4), intent(in)  :: iarray1
        integer(i4), intent(out) :: iarray2
     end subroutine all_reduce_sum_i41

     subroutine all_reduce_sum_r8(farray1,farray2,ldim)
        use my_kinddefs
        implicit none
        real(r8)   , intent(in)  :: farray1(ldim)
        real(r8)   , intent(out) :: farray2(ldim)
        integer(i4), intent(in)  :: ldim
     end subroutine all_reduce_sum_r8

     subroutine all_reduce_sum_r81(farray1,farray2)
        use my_kinddefs
        implicit none
        real(r8), intent(in)  :: farray1
        real(r8), intent(out) :: farray2
     end subroutine all_reduce_sum_r81

  end interface all_reduce_sum

end module all_reduce_sum_mod

!-------------------------------------------------------------------------------
subroutine all_reduce_sum_i4(iarray1,iarray2,ldim)
!Version for array of i4 integers (ldim of them) ...

      use my_kinddefs
      implicit none
#include "mympif.h"

      integer(i4), intent(in)  :: iarray1(ldim)
      integer(i4), intent(out) :: iarray2(ldim)
      integer(i4), intent(in)  :: ldim
!--Tmps
      integer(i4) ::   i,ierr

#ifdef MPI_ON
      call MPI_ALLREDUCE(iarray1,iarray2,ldim,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      do i=1,ldim
       iarray2(i)  = iarray1(i) 
      enddo
#endif

end subroutine all_reduce_sum_i4
!-------------------------------------------------------------------------------
subroutine all_reduce_sum_i41(iarray1,iarray2)
!Version for single i4 integer...

      use my_kinddefs
      implicit none
#include "mympif.h"

      integer(i4), intent(in)  :: iarray1
      integer(i4), intent(out) :: iarray2
!--Tmps
      integer(i4) ::   ldim
      integer(i4) ::   i,ierr

#ifdef MPI_ON
      ldim     = 1
      call MPI_ALLREDUCE(iarray1,iarray2,ldim,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      iarray2  = iarray1 
#endif

end subroutine all_reduce_sum_i41
!-------------------------------------------------------------------------------
subroutine all_reduce_sum_r8(farray1,farray2,ldim)
!Version for array of r8 reals (ldim of them) ...

      use my_kinddefs
      implicit none
#include "mympif.h"

      real(r8)   , intent(in)  :: farray1(ldim)
      real(r8)   , intent(out) :: farray2(ldim)
      integer(i4), intent(in)  :: ldim
!--Tmps
      integer(i4) ::   i,ierr

#ifdef MPI_ON
      call MPI_ALLREDUCE(farray1,farray2,ldim,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      do i=1,ldim
       farray2(i)  = farray1(i)
      enddo
#endif

end subroutine all_reduce_sum_r8
!-------------------------------------------------------------------------------
subroutine all_reduce_sum_r81(farray1,farray2)
!Version for single r8 real...

      use my_kinddefs
      implicit none
#include "mympif.h"

      real(r8)   , intent(in)  :: farray1
      real(r8)   , intent(out) :: farray2
!--Tmps
      integer(i4) ::   i,ierr
      integer(i4) ::   ldim

#ifdef MPI_ON
      ldim        = 1
      call MPI_ALLREDUCE(farray1,farray2,ldim,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
       farray2  = farray1
#endif

end subroutine all_reduce_sum_r81
!-------------------------------------------------------------------------------
