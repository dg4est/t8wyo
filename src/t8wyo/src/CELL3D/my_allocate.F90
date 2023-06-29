module my_allocate_mod
!*****************************************************************************80
!> \brief my_allocate: Allocate and perform other functions for all pointer types
!!
!> \ details
!! This module contains all dynamic memory allocation for pointers.
!! The following steps are taken:
!! 1: Pointer is checked for association.  If it is associated, error is flagged. Otherwise execution continues.
!!    This feature can be disabled by setting POINTER_ASSOCIATION_CHECK= 0 in io_params.f90
!! 2: Size of allocation is checked.  If it is < 0, error is flagged.
!!    If it is = 0, size is reset internally = 1 since size 0 allocation does not result in pointer association.
!! 3: Pointer is allocated (and thus associated)
!! 4: Pointer/array values are initialize = 0
!! 5: Message is written to log and procedure exits.
!!    Note: character allocates currently require beginning and ending sizes and option to do this for i4 is included.
!*****************************************************************************80

  interface my_allocate

     subroutine my_allocate_i4(iarray,isize,ctag)
        use my_kinddefs
        integer(i4),pointer :: iarray(:)
        integer(i4)         ::  isize
        character(*)        ::  ctag
     end subroutine my_allocate_i4

     subroutine my_allocate_i04(iarray,isize0,isize,ctag)
        use my_kinddefs
        integer(i4),pointer :: iarray(:)
        integer(i4)         ::  isize0
        integer(i4)         ::  isize
        character(*)        ::  ctag
     end subroutine my_allocate_i04

     subroutine my_allocate_i8(iarray,isize,ctag)
        use my_kinddefs
        integer(i8),pointer :: iarray(:)
        integer(i4)         ::  isize
        character(*)        ::  ctag
     end subroutine my_allocate_i8

     subroutine my_allocate_r4(farray,isize,ctag)
        use my_kinddefs
        real(r4), pointer   :: farray(:)
        integer(i4)         :: isize
        character(*)        ::  ctag
     end subroutine my_allocate_r4

     subroutine my_allocate_r8(farray,isize,ctag)
        use my_kinddefs
        real(r8), pointer   :: farray(:)
        integer(i4)         :: isize
        character(*)        ::  ctag
     end subroutine my_allocate_r8

     subroutine my_allocate_ch20(carray,isize0,isize,ctag)
        use my_kinddefs
        character*20,pointer:: carray(:)
        integer(i4)         ::  isize0
        integer(i4)         ::  isize
        character(*)        ::  ctag
     end subroutine my_allocate_ch20

     subroutine my_allocate_int_ptr_array(iptr,isize,ctag)
        use my_kinddefs
        use my_typedefs
        type(int_ptr_array),pointer ::  iptr(:)
        integer(i4)                 ::  isize
        character(*)                ::  ctag
     end subroutine my_allocate_int_ptr_array

     subroutine my_allocate_real_ptr_array(fptr,isize,ctag)
        use my_kinddefs
        use my_typedefs
        type(real_ptr_array),pointer ::  fptr(:)
        integer(i4)                 ::  isize
        character(*)                ::  ctag
     end subroutine my_allocate_real_ptr_array


     subroutine my_allocate_mpi_array(mpi,isize,ctag)
        use my_kinddefs
        use my_typedefs
        type(mpi_array),pointer     ::  mpi(:)
        integer(i4)                 ::  isize
        character(*)                ::  ctag
     end subroutine my_allocate_mpi_array



  end interface my_allocate
end module my_allocate_mod

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_allocate_i4(iarray,isize_in,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  integer(i4), pointer              :: iarray(:)   !(intentout)!
  integer(i4),          INTENT(IN)  :: isize_in
  character(*),         INTENT(IN)  :: ctag

  integer(i4) :: isize
  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_ASSOCIATION_CHECK==1 .and. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif


  if (isize_in < 0) then
    write(iwrit ,697) ctag                                  !Failed
    call abort_all
  endif


  isize = max(isize_in,1)                                   !Allocate minimum size=1
                                                            !Otherwise, pointers of size=0
                                                            !cannot be passed to subroutines
                                                      
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

  allocate(iarray(isize), STAT = istatus)                   !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
    iarray = 0                                              !Optionally, initialize to 0
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
697 format( '(*)ERROR: NEGATIVE SIZE ALLOCATION   FOR : ',a)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_allocate_i04(iarray,isize0,isize,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  integer(i4), pointer              :: iarray(:)   !(intent(out)!
  integer(i4),          INTENT(IN)  :: isize0
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (isize < isize0) then 
    write(iwrit ,697) ctag                                  !Failed
    call abort_all
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

  allocate(iarray(isize0:isize), STAT = istatus)            !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
  iarray(isize0:isize) = 0                                  !Optionally initialize to 0
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
697 format( '(*)ERROR: NEGATIVE SIZE ALLOCATION   FOR : ',a)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_i04
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_allocate_i8(iarray,isize_in,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  integer(i8), pointer              :: iarray(:)  !(intent(out)!
  integer(i4),          INTENT(IN)  :: isize_in
  character(*),         INTENT(IN)  :: ctag

  integer(i4) :: isize
  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_ASSOCIATION_CHECK==1 .and. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif

  if (isize_in < 0) then
    write(iwrit ,697) ctag                                  !Failed
    call abort_all
  endif


  isize = max(isize_in,1)                                   !Allocate minimum size=1
                                                            !Otherwise, pointers of size=0

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

  allocate(iarray(isize), STAT = istatus)                   !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
   iarray = 0                                               !Optionally initialize to 0
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
697 format( '(*)ERROR: NEGATIVE SIZE ALLOCATION   FOR : ',a)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_i8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_allocate_r4(farray,isize_in,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  real(r4)   , pointer              :: farray(:)  !(intent(out)!
  integer(i4),          INTENT(IN)  :: isize_in
  character(*),         INTENT(IN)  :: ctag

  integer(i4) :: isize
  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_ASSOCIATION_CHECK==1 .and. associated(farray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif

  if (isize_in < 0) then
    write(iwrit ,697) ctag                                  !Failed
    call abort_all
  endif


  isize = max(isize_in,1)                                   !Allocate minimum size=1
                                                            !Otherwise, pointers of size=0
                                                            !cannot be passed to subroutines

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

  allocate(farray(isize), STAT = istatus)                   !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
   farray(isize) = 0.0                                      !Optionally initialize to 0
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
697 format( '(*)ERROR: NEGATIVE SIZE ALLOCATION   FOR : ',a)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_allocate_r8(farray,isize_in,ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  real(r8)   , pointer              :: farray(:)   !intent(out)!
  integer(i4),          INTENT(IN)  :: isize_in
  character(*),         INTENT(IN)  :: ctag

  integer(i4) :: isize
  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units

  if (POINTER_ASSOCIATION_CHECK==1 .and. associated(farray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif

  if (isize_in < 0) then
    write(iwrit ,697) ctag                                  !Failed
    call abort_all
  endif


  isize = max(isize_in,1)                                   !Allocate minimum size=1
                                                            !Otherwise, pointers of size=0
                                                            !cannot be passed to subroutines
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

  allocate(farray(isize), STAT = istatus)                   !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
   farray = 0.0                                             !Optionally initialize to 0.0
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
697 format( '(*)ERROR: NEGATIVE SIZE ALLOCATION   FOR : ',a)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_allocate_ch20(carray,isize0,isize,ctag)

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

  if (POINTER_ASSOCIATION_CHECK==1 .and. associated(carray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif

  if (isize < isize0) then
    write(iwrit ,697) ctag                                  !Failed
    call abort_all
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

  allocate(carray(isize0:isize), STAT = istatus)                   !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
! carray(isize0:isize) = ' '                                !Optionally initialize to ' '
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
697 format( '(*)ERROR: NEGATIVE SIZE ALLOCATION   FOR : ',a)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_ch20
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine my_allocate_int_ptr_array(iarray,isize,ctag)

  use my_kinddefs
  use my_typedefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  type(int_ptr_array), pointer      :: iarray(:)   !(intentout)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units
  integer(i4) :: i

  if (POINTER_ASSOCIATION_CHECK==1 .and. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif

  if (isize .le. 0) then
    allocate(iarray(1)    , STAT = istatus)                   !Allocate minimum size=1
    if (istatus .eq. 0) then                                  !Otherwise, pointers of size=0
      return                                                  !cannot be passed to subroutines
    else
!     write(io_mem,699) ctag                                  !Failed (OMIT because other PE may overwrite entire hist file)
      write(iwrit ,699) ctag                                  !Failed
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

  allocate(iarray(isize), STAT = istatus)                   !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
   do i=1,isize                                             !Must nullify ptrs
     nullify(iarray(i)%ptr)                                 !in this defined type
   enddo                                                    !to ensure they can be allocated later
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_int_ptr_array
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_allocate_real_ptr_array(iarray,isize,ctag)

  use my_kinddefs
  use my_typedefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"
!
  type(real_ptr_array), pointer     :: iarray(:)   !(intentout)!
  integer(i4),          INTENT(IN)  :: isize
  character(*),         INTENT(IN)  :: ctag

  integer(i8) :: isize_mem
  integer(i4) :: istatus
  real(r8)    :: fsize_mem
  character*7 :: mem_units
  integer(i4) :: i

  if (POINTER_ASSOCIATION_CHECK==1 .and. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif

  if (isize .le. 0) then
    allocate(iarray(1)    , STAT = istatus)                   !Allocate minimum size=1
    if (istatus .eq. 0) then                                  !Otherwise, pointers of size=0
      return                                                  !cannot be passed to subroutines
    else
!     write(io_mem,699) ctag                                  !Failed (OMIT because other PE may overwrite entire hist file)
      write(iwrit ,699) ctag                                  !Failed
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

  allocate(iarray(isize), STAT = istatus)                   !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
   do i=1,isize                                             !Must nullify ptrs
     nullify(iarray(i)%ptr)                                 !in this defined type
   enddo                                                    !to ensure they can be allocated later
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_real_ptr_array
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine my_allocate_mpi_array(iarray,isize,ctag)

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
  integer(i4) :: i

  if (POINTER_ASSOCIATION_CHECK==1 .and. associated(iarray)) then
    write(io_mem,698) ctag
    write(iwrit,698) ctag
    call abort_all
  endif

  if (isize .le. 0) then
    allocate(iarray(1)    , STAT = istatus)                   !Allocate minimum size=1
    if (istatus .eq. 0) then                                  !Otherwise, pointers of size=0
      return                                                  !cannot be passed to subroutines
    else
!     write(io_mem,699) ctag                                  !Failed (OMIT because other PE may overwrite entire hist file)
      write(iwrit ,699) ctag                                  !Failed
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

  allocate(iarray(isize), STAT = istatus)                   !Perform actual allocate...

  if (istatus .eq. 0) then
    ntotal_dyn_mem = ntotal_dyn_mem + isize_mem             !Sucessful
  else
!   write(io_mem,699) ctag                                  !Failed
    write(iwrit ,699) ctag                                  !Failed
    call abort_all
  endif

!----------------------------------------------------------------------
   do i=1,isize                            !Nullify all ptrs in this defined type
      nullify(iarray(i)%ptr%iproc_send)       !to ensure these can be used later
      nullify(iarray(i)%ptr%iproc_recv)       !in my_allocate
      nullify(iarray(i)%ptr%ipntr_send)
      nullify(iarray(i)%ptr%ipntr_recv)
      nullify(iarray(i)%ptr%ibuff_send)
      nullify(iarray(i)%ptr%ibuff_recv)
      nullify(iarray(i)%ptr%fbuff_send)
      nullify(iarray(i)%ptr%fbuff_recv)
      nullify(iarray(i)%ptr%ilocal_send)
      nullify(iarray(i)%ptr%ilocal_recv)
      nullify(iarray(i)%ptr%msgid)
      nullify(iarray(i)%ptr%istatus)
   enddo
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

601 format( '(+)Allocating memory for   : ',a)
602 format( '(+)Requested Memory to be allocated  :            ',e15.4,a)
603 format( '(+)Memory Allocated Successfully  : TOTAL MEMORY: ',e15.4,a,/)
698 format( '(*)ERROR: POINTER ALREADY ASSOCIATED FOR : ',a)
699 format( '(*)MEMORY ALLOCATION FAILED FOR  : ',a)

end subroutine my_allocate_mpi_array
!------------------------------------------------------------------------------
