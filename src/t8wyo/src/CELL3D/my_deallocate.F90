module my_deallocate_mod

!***********************************************
!***********************************************
!--Should probably nullify all ptrs on deallocation...*** (5/25/12)
!--This has been done. (2/11/13)
!--Another issue: Check all member in arrays of pts to see if they
!--are null.  If not, this may correspond to mem that is no longer accessible.
!!-Perhaps may this an option in my_deallocate (since deallocating parent may be desired in some cases).
!***********************************************
!***********************************************
  interface my_deallocate

     subroutine my_deallocate_i4(iarray,isize,ctag)
        use my_kinddefs
        integer(i4),pointer :: iarray(:)
        integer(i4)         ::  isize
        character(*)        ::  ctag
     end subroutine my_deallocate_i4

     subroutine my_deallocate_i04(iarray,isize0,isize,ctag)
        use my_kinddefs
        integer(i4),pointer :: iarray(:)
        integer(i4)         ::  isize0
        integer(i4)         ::  isize
        character(*)        ::  ctag
     end subroutine my_deallocate_i04

     subroutine my_deallocate_i8(iarray,isize,ctag)
        use my_kinddefs
        integer(i8),pointer :: iarray(:)
        integer(i4)         ::  isize
        character(*)        ::  ctag
     end subroutine my_deallocate_i8

     subroutine my_deallocate_r4(farray,isize,ctag)
        use my_kinddefs
        real(r4), pointer   :: farray(:)
        integer(i4)         :: isize
        character(*)        ::  ctag
     end subroutine my_deallocate_r4

     subroutine my_deallocate_r8(farray,isize,ctag)
        use my_kinddefs
        real(r8), pointer   :: farray(:)
        integer(i4)         :: isize
        character(*)        ::  ctag
     end subroutine my_deallocate_r8

     subroutine my_deallocate_ch20(carray,isize0,isize,ctag)
        use my_kinddefs
        character*20,pointer:: carray(:)
        integer(i4)         ::  isize0
        integer(i4)         ::  isize
        character(*)        ::  ctag
     end subroutine my_deallocate_ch20

     subroutine my_deallocate_int_ptr_array(iptr,isize,ctag)
        use my_kinddefs
        use my_typedefs
        type(int_ptr_array),pointer ::  iptr(:)
        integer(i4)                 ::  isize
        character(*)                ::  ctag
     end subroutine my_deallocate_int_ptr_array

     subroutine my_deallocate_real_ptr_array(fptr,isize,ctag)
        use my_kinddefs
        use my_typedefs
        type(real_ptr_array),pointer::  fptr(:)
        integer(i4)                 ::  isize
        character(*)                ::  ctag
     end subroutine my_deallocate_real_ptr_array


     subroutine my_deallocate_mpi_array(mpi,isize,ctag)
        use my_kinddefs
        use my_typedefs
        type(mpi_array),pointer     ::  mpi(:)
        integer(i4)                 ::  isize
        character(*)                ::  ctag
     end subroutine my_deallocate_mpi_array



  end interface my_deallocate
end module my_deallocate_mod

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_i4(iarray,isize,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  integer(i4), pointer              :: iarray(:)   !intent(out)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_DEASSOCIATION_CHECK == 1) then
  if (.not. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif
  endif

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(iarray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(iarray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed : OMIT since this could be another PE overwriting file
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif

  isize_mem = isize * i4_bytesize

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(iarray, STAT = istatus)                   !Perform actual deallocate...
    nullify(iarray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_i04(iarray,isize0,isize,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  integer(i4), pointer              :: iarray(:)  !intent(out)!
  integer(i4),          INTENT(IN)  :: isize0
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_DEASSOCIATION_CHECK == 1) then
  if (.not. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif
  endif

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(iarray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(iarray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif


  isize_mem = (isize - isize0 + 1) * i4_bytesize

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(iarray, STAT = istatus)            !Perform actual deallocate...
  nullify(iarray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_i04
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_i8(iarray,isize,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  integer(i8), pointer              :: iarray(:)   !intent(out)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_DEASSOCIATION_CHECK == 1) then
  if (.not. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif
  endif

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(iarray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(iarray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif


  isize_mem = isize * i8_bytesize

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(iarray, STAT = istatus)                   !Perform actual deallocate...
  nullify(iarray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_i8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_r4(farray,isize,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  real(r4)   , pointer              :: farray(:)   !intent(out)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_DEASSOCIATION_CHECK == 1) then
  if (.not. associated(farray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif
  endif

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(farray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(farray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif


  isize_mem = isize * r4_bytesize

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(farray, STAT = istatus)                   !Perform actual deallocate...
  nullify(farray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_r8(farray,isize,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  real(r8)   , pointer              :: farray(:)   !intent(out)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_DEASSOCIATION_CHECK == 1) then
  if (.not. associated(farray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif
  endif

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(farray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(farray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif


  isize_mem = isize * r8_bytesize

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(farray, STAT = istatus)                   !Perform actual deallocate...
  nullify(farray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_ch20(carray,isize0,isize,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  character*20,pointer              :: carray(:)   !intent(out)!
  integer(i4),          INTENT(IN)  :: isize0
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(carray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(carray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif


  isize_mem = (isize - isize0 + 1) * i4_bytesize*20

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(carray, STAT = istatus)                   !Perform actual deallocate...
  nullify(carray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_ch20
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_int_ptr_array(iarray,isize,ctag)

  use my_kinddefs
  use my_typedefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  type(int_ptr_array),pointer       :: iarray(:)   !intent(out)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_DEASSOCIATION_CHECK == 1) then
  if (.not. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif
  endif

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(iarray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(iarray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed : OMIT since this could be another PE overwriting file
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif

  isize_mem = isize * int_ptr_size  !<--Best estimate of size of ptrs in this defined type

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(iarray, STAT = istatus)                   !Perform actual deallocate...
  nullify(iarray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_int_ptr_array
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_real_ptr_array(iarray,isize,ctag)

  use my_kinddefs
  use my_typedefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  type(real_ptr_array),pointer      :: iarray(:)   !intent(out)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_DEASSOCIATION_CHECK == 1) then
  if (.not. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif
  endif

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(iarray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(iarray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed : OMIT since this could be another PE overwriting file
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif

  isize_mem = isize * ireal_ptr_size  !<--Best estimate of size of ptrs in this defined type

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(iarray, STAT = istatus)                   !Perform actual deallocate...
  nullify(iarray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_real_ptr_array
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_deallocate_mpi_array(iarray,isize,ctag)

  use my_kinddefs
  use my_typedefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  type(mpi_array), pointer          :: iarray(:)   !(intentout)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_DEASSOCIATION_CHECK == 1) then
  if (.not. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif
  endif

  if (isize .le. 0) then                                 !For zero size array...
    deallocate(iarray, STAT = istatus)                   !deallocate (was allocated as size=1)
    nullify(iarray)
    if (istatus .eq. 0) then
      return
    else
!     write(io_mem,699) ctag                             !Failed : OMIT since this could be another PE overwriting file
      write(iwrit ,699) ctag                             !Failed
      call abort_all
    endif
  endif

  isize_mem = isize * mpi_ptr_size  !<--Best estimate of size of ptrs in this defined type

  if     (isize_mem .lt. 1000000) then
     fsize_mem = dble(isize_mem)
     mem_units ='  Bytes'
  elseif (isize_mem .lt. 1000000000) then
     fsize_mem = dble(isize_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(isize_mem)/1000000000.
     mem_units =' GBytes'
  endif

  if (id_proc .eq. io_mem_proc) then
   write(io_mem,601) ctag
   write(io_mem,602) fsize_mem,mem_units
   call flush(io_mem)
  endif


!----------------------------------------------------------------------

  deallocate(iarray, STAT = istatus)                   !Perform actual deallocate...
  nullify(iarray)

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem - isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------

  if     (ntotal_dyn_mem .lt. 1000000) then
     fsize_mem = dble(ntotal_dyn_mem)
     mem_units ='  Bytes'
  elseif (ntotal_dyn_mem .lt. 1000000000) then
     fsize_mem = dble(ntotal_dyn_mem)/1000000.
     mem_units =' MBytes'
  else
     fsize_mem = dble(ntotal_dyn_mem)/1000000000.
     mem_units =' GBytes'
  endif


  if (id_proc .eq. io_mem_proc) then
   write(io_mem,603) fsize_mem,mem_units
   call flush(io_mem)
  endif

601 format('(-)Deallocating memory for : ',a)
602 format('(-)Requested Memory to be deallocated:            ',e15.4,a)
603 format('(-)Memory Deallocated Successfully: TOTAL MEMORY: ',e15.4,a,/)
698 format('(*)ERROR: POINTER NOT ASSOCIATED BEFORE DEALLOCATION FOR : ',a)
699 format('(*)MEMORY DEALLOCATION FAILED FOR: ',a)

end subroutine my_deallocate_mpi_array
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
