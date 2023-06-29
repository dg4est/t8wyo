!-------------------------------------------------------------------------------
module mpi_to_gpt_mod

!--Copy real point values to gpt locations
!--Allows for different input and output arrays and dimensions
!--Arrays may be either: 1D pointers or 2D arrays (i4 and i42 respectively)
!--Most general case when real and gpt buffers are different arrays  (i42 and i420 routines)
!--Real points: iarray(:)     or farray(lblock,nnode)
!--Gpt  points: iarray_new(:) or farray_new(lblock_new,nnode_new)
!--Most common case will be when these are same arrays 
!--and simpler calling sequence is provided (i40 and i4 routines)
!--Thus, total of 4 calling sequences:
!--1: single array,  dimensioned as 1D pointer
!--2: single array,  dimensioned as 2D array
!--3: two    arrays, dimensioned as 1D pointers
!--4: two    arrays, dimensioned as 2D arrays  <--Detailed implementation corresponds to this one (i42)
!--Repeat for r8 routines

  interface mpi_to_gpt

!--Integer version with iarray ONLY as 1D pointers...
    subroutine mpi_to_gpt_i400(nnode,iarray,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      integer(i4),          intent(inout) :: iarray(nnode)
      type(mpi_schedule),   intent(inout) :: mpi0
    end subroutine mpi_to_gpt_i400

!--Integer version with iarray ONLY as 1D pointers...
    subroutine mpi_to_gpt_i40(nnode,lblock,kblock,iarray,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      integer(i4),          intent(in)    :: lblock
      integer(i4),          intent(in)    :: kblock
      integer(i4), pointer                :: iarray(:)  !intent(inout)!
      type(mpi_schedule),   intent(inout) :: mpi0
    end subroutine mpi_to_gpt_i40

!--Integer version with iarray ONLY as 2D array...
    subroutine mpi_to_gpt_i4(nnode,lblock,kblock,iarray,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      integer(i4),          intent(in)    :: lblock
      integer(i4),          intent(in)    :: kblock
      integer(i4),          intent(inout) :: iarray(lblock,nnode)
      type(mpi_schedule),   intent(inout) :: mpi0
    end subroutine mpi_to_gpt_i4


!--Integer version with iarray and iarray_new as 1D pointers...
    subroutine mpi_to_gpt_i420(nnode,lblock,kblock,iarray,                    &
                               nnode_new,lblock_new,kblock_new,iarray_new,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      integer(i4),          intent(in)    :: lblock
      integer(i4),          intent(in)    :: kblock
      integer(i4), pointer                :: iarray(:)  !intent(in)!
      integer(i4),          intent(in)    :: nnode_new
      integer(i4),          intent(in)    :: lblock_new
      integer(i4),          intent(in)    :: kblock_new
      integer(i4), pointer                :: iarray_new(:) !intent(inout)!
      type(mpi_schedule),   intent(inout) :: mpi0
    end subroutine mpi_to_gpt_i420


!--Integer version with iarray and iarray_new as 2D arrays...
    subroutine mpi_to_gpt_i42(nnode,lblock,kblock,iarray,                    &
                              nnode_new,lblock_new,kblock_new,iarray_new,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4), intent(in)    :: nnode
      integer(i4), intent(in)    :: lblock
      integer(i4), intent(in)    :: kblock
      integer(i4), intent(inout) :: iarray(lblock,nnode)
      integer(i4), intent(in)    :: nnode_new
      integer(i4), intent(in)    :: lblock_new
      integer(i4), intent(in)    :: kblock_new
      integer(i4), intent(inout) :: iarray_new(lblock_new,nnode_new)
      type(mpi_schedule), intent(inout) :: mpi0
    end subroutine mpi_to_gpt_i42

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------

!--REAL8   version with iarray ONLY as 1D pointers...
    subroutine mpi_to_gpt_r800(nnode,iarray,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      real(r8),             intent(inout) :: iarray(nnode)
      type(mpi_schedule),   intent(inout) :: mpi0
    end subroutine mpi_to_gpt_r800

!--REAL8   version with iarray ONLY as 1D pointers...
    subroutine mpi_to_gpt_r80(nnode,lblock,kblock,iarray,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      integer(i4),          intent(in)    :: lblock
      integer(i4),          intent(in)    :: kblock
      real(r8),    pointer                :: iarray(:) !intent(inout)!
      type(mpi_schedule),   intent(inout) :: mpi0
    end subroutine mpi_to_gpt_r80

!--REAL8   version with iarray ONLY as 2D array...
    subroutine mpi_to_gpt_r8(nnode,lblock,kblock,iarray,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      integer(i4),          intent(in)    :: lblock
      integer(i4),          intent(in)    :: kblock
      real(r8),             intent(inout) :: iarray(lblock,nnode)
      type(mpi_schedule),   intent(inout) :: mpi0
    end subroutine mpi_to_gpt_r8


!--REAL8   version with iarray and iarray_new as 1D pointers...
    subroutine mpi_to_gpt_r820(nnode,lblock,kblock,iarray,                    &
                               nnode_new,lblock_new,kblock_new,iarray_new,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      integer(i4),          intent(in)    :: lblock
      integer(i4),          intent(in)    :: kblock
      real(r8),    pointer                :: iarray(:) !intent(in)!
      integer(i4),          intent(in)    :: nnode_new
      integer(i4),          intent(in)    :: lblock_new
      integer(i4),          intent(in)    :: kblock_new
      real(r8),    pointer                :: iarray_new(:) !intent(inout)!
      type(mpi_schedule),   intent(inout) :: mpi0
    end subroutine mpi_to_gpt_r820


!--REAL8   version with iarray and iarray_new as 2D arrays...
    subroutine mpi_to_gpt_r82(nnode,lblock,kblock,iarray,                    &
                              nnode_new,lblock_new,kblock_new,iarray_new,mpi0)
      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4), intent(in)    :: nnode
      integer(i4), intent(in)    :: lblock
      integer(i4), intent(in)    :: kblock
      real(r8),    intent(inout) :: iarray(lblock,nnode)
      integer(i4), intent(in)    :: nnode_new
      integer(i4), intent(in)    :: lblock_new
      integer(i4), intent(in)    :: kblock_new
      real(r8),    intent(inout) :: iarray_new(lblock_new,nnode_new)
      type(mpi_schedule), intent(inout) :: mpi0
    end subroutine mpi_to_gpt_r82

  end interface mpi_to_gpt

end module mpi_to_gpt_mod
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
 subroutine mpi_to_gpt_i400(nnode,iarray,mpi0)

      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),        intent(in)      :: nnode
      integer(i4),        intent(inout)   :: iarray(nnode)
      type(mpi_schedule), intent(inout)   :: mpi0
      integer(i4)  :: lblock

      lblock = 1
      call mpi_to_gpt_i42(nnode,lblock,lblock,iarray,    &
                          nnode,lblock,lblock,iarray,mpi0)
end subroutine mpi_to_gpt_i400
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
 subroutine mpi_to_gpt_i40(nnode,lblock,kblock,iarray,mpi0)

      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4), intent(in)             :: nnode
      integer(i4), intent(in)             :: lblock
      integer(i4), intent(in)             :: kblock
      integer(i4), pointer                :: iarray(:) !intent(inout)!
      type(mpi_schedule), intent(inout)   :: mpi0

      call mpi_to_gpt_i42(nnode,lblock,kblock,iarray,    &
                          nnode,lblock,kblock,iarray,mpi0)
end subroutine mpi_to_gpt_i40
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 subroutine mpi_to_gpt_i4(nnode,lblock,kblock,iarray,mpi0)

      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4), intent(in)           :: nnode
      integer(i4), intent(in)           :: lblock
      integer(i4), intent(in)           :: kblock
      integer(i4), intent(inout)        :: iarray(lblock,nnode)
      type(mpi_schedule), intent(inout) :: mpi0

      call mpi_to_gpt_i42(nnode,lblock,kblock,iarray,    &
                          nnode,lblock,kblock,iarray,mpi0)
end subroutine mpi_to_gpt_i4
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!--Integer version with iarray and iarray_new as 1D pointers...
 subroutine mpi_to_gpt_i420(nnode,lblock,kblock,iarray,                    &
                            nnode_new,lblock_new,kblock_new,iarray_new,mpi0)

      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4), intent(in)             :: nnode
      integer(i4), intent(in)             :: lblock
      integer(i4), intent(in)             :: kblock
      integer(i4), pointer                :: iarray(:)  !intent(in)!
      integer(i4), intent(in)             :: nnode_new
      integer(i4), intent(in)             :: lblock_new
      integer(i4), intent(in)             :: kblock_new
      integer(i4), pointer                :: iarray_new(:)  !intent(inout)!
      type(mpi_schedule), intent(inout)   :: mpi0

      call mpi_to_gpt_i42(nnode,lblock,kblock,iarray,                    &
                          nnode_new,lblock_new,kblock_new,iarray_new,mpi0)
end subroutine mpi_to_gpt_i420

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!--Integer version with iarray and iarray_new as 2D arrays...
  subroutine mpi_to_gpt_i42(nnode,lblock,kblock,iarray,                    &
                            nnode_new,lblock_new,kblock_new,iarray_new,mpi0)

!--Reverse exchange:  Send from Real Pts to Gpts
!--MPI Sends recv buffers and receives send buffers...


      use my_kinddefs
      use mpi_schedule_typedef
      use mp_stuff
      use my_allocate_mod
      use my_deallocate_mod
      use my_mpi_barrier_mod
      implicit none
#include "mympif.h"

      integer(i4), intent(in)           :: nnode
      integer(i4), intent(in)           :: lblock
      integer(i4), intent(in)           :: kblock
      integer(i4), intent(inout)        :: iarray(lblock,nnode)
      integer(i4), intent(in)           :: nnode_new
      integer(i4), intent(in)           :: lblock_new
      integer(i4), intent(in)           :: kblock_new
      integer(i4), intent(inout)        :: iarray_new(lblock_new,nnode_new)
      type(mpi_schedule), intent(inout) :: mpi0
!--Tmps
      integer(i4) :: ii,in,iproc,ip,iip,iip2
      integer(i4) :: kptr,ncount


!-----------------------------------------------------------------
!--Allocate buffer arrays for mpi communication
       call my_allocate(mpi0%ibuff_send,mpi0%nbuff_send,'mpi_to_gpt_i4:mpi0%ibuff_send')
       call my_allocate(mpi0%ibuff_recv,mpi0%nbuff_recv,'mpi_to_gpt_i4:mpi0%ibuff_recv')
!-----------------------------------------------------------------
!------------------------------------------
!--Fill send buffers
      do ii=1,mpi0%nbuff_recv
         mpi0%ibuff_recv(ii) = iarray(kblock,mpi0%ilocal_recv(ii))
      enddo
!------------------------------------------
      mpi0%itype_recv = 10000
      mpi0%itype_send = MPI_ANY_TAG

!--Post RECVs
      do iip=1,mpi0%nproc_send
        iproc  = mpi0%iproc_send(iip)-1
        kptr   = mpi0%ipntr_send(iip)
        ncount =  mpi0%ipntr_send(iip+1) - mpi0%ipntr_send(iip)
        call MPI_IRECV(mpi0%ibuff_send(kptr),ncount,MPI_INTEGER,   &
                       iproc,mpi0%itype_send,MPI_COMM_WORLD,       &
                       mpi0%msgid(iip),mpi0%ierr_send)

      enddo

!--Perform SENDS
      do iip=1,mpi0%nproc_recv
        iip2   = iip + mpi0%nproc_send
        iproc  = mpi0%iproc_recv(iip)-1
        kptr   = mpi0%ipntr_recv(iip)
        ncount =  mpi0%ipntr_recv(iip+1) - mpi0%ipntr_recv(iip)
        call MPI_ISEND(mpi0%ibuff_recv(kptr),ncount,MPI_INTEGER,  &
                       iproc,mpi0%itype_recv,MPI_COMM_WORLD,      &
                       mpi0%msgid(iip2),mpi0%ierr_recv)

      enddo

!--Wait for All Messages to Arrive
      call MPI_WAITALL(mpi0%nproc_all,mpi0%msgid,mpi0%istatus,mpi0%ierr)


!--Empty Recv buffers
      do ii=1,mpi0%nbuff_send
        iarray_new(kblock_new,mpi0%ilocal_send(ii)) = mpi0%ibuff_send(ii)
      enddo

      call my_mpi_barrier(1) !Use barrier after waitall

!-----------------------------------------------------------------
!--Deallocate buffer arrays for mpi communication
       call my_deallocate(mpi0%ibuff_send,mpi0%nbuff_send,'mpi_to_gpt_i4:mpi0%ibuff_send')
       call my_deallocate(mpi0%ibuff_recv,mpi0%nbuff_recv,'mpi_to_gpt_i4:mpi0%ibuff_recv')
!-----------------------------------------------------------------

end subroutine mpi_to_gpt_i42




!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
 subroutine mpi_to_gpt_r800(nnode,iarray,mpi0)

      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4),          intent(in)    :: nnode
      real(r8),             intent(inout) :: iarray(nnode)
      type(mpi_schedule),   intent(inout) :: mpi0
      integer(i4) :: lblock

      lblock = 1
      call mpi_to_gpt_r82(nnode,lblock,lblock,iarray,    &
                          nnode,lblock,lblock,iarray,mpi0)
end subroutine mpi_to_gpt_r800
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
 subroutine mpi_to_gpt_r80(nnode,lblock,kblock,iarray,mpi0)

      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4), intent(in)             :: nnode
      integer(i4), intent(in)             :: lblock
      integer(i4), intent(in)             :: kblock
      real(r8),    pointer                :: iarray(:)  !intent(inout)!
      type(mpi_schedule), intent(inout)   :: mpi0

      call mpi_to_gpt_r82(nnode,lblock,kblock,iarray,    &
                          nnode,lblock,kblock,iarray,mpi0)
end subroutine mpi_to_gpt_r80
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 subroutine mpi_to_gpt_r8(nnode,lblock,kblock,iarray,mpi0)

      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4), intent(in)           :: nnode
      integer(i4), intent(in)           :: lblock
      integer(i4), intent(in)           :: kblock
      real(r8),    intent(inout)        :: iarray(lblock,nnode)
      type(mpi_schedule), intent(inout) :: mpi0

      call mpi_to_gpt_r82(nnode,lblock,kblock,iarray,    &
                          nnode,lblock,kblock,iarray,mpi0)
end subroutine mpi_to_gpt_r8
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!--Integer version with iarray and iarray_new as 1D pointers...
 subroutine mpi_to_gpt_r820(nnode,lblock,kblock,iarray,                    &
                            nnode_new,lblock_new,kblock_new,iarray_new,mpi0)

      use my_kinddefs
      use mpi_schedule_typedef
      implicit none
      integer(i4), intent(in)             :: nnode
      integer(i4), intent(in)             :: lblock
      integer(i4), intent(in)             :: kblock
      real(r8),    pointer                :: iarray(:)  !intent(in)!
      integer(i4), intent(in)             :: nnode_new
      integer(i4), intent(in)             :: lblock_new
      integer(i4), intent(in)             :: kblock_new
      real(r8),    pointer                :: iarray_new(:) !intent(inout)!
      type(mpi_schedule), intent(inout)   :: mpi0

      call mpi_to_gpt_r82(nnode,lblock,kblock,iarray,                    &
                          nnode_new,lblock_new,kblock_new,iarray_new,mpi0)
end subroutine mpi_to_gpt_r820

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine mpi_to_gpt_r82(nnode,lblock,kblock,farray,                    &
                            nnode_new,lblock_new,kblock_new,farray_new,mpi0)

!--Reverse exchange:  Send from Real Pts to Gpts
!--MPI Sends recv buffers and receives send buffers...


      use my_kinddefs
      use mpi_schedule_typedef
      use mp_stuff
      use my_allocate_mod
      use my_deallocate_mod
      use my_mpi_barrier_mod
      implicit none
#include "mympif.h"

      integer(i4), intent(in)    :: nnode
      integer(i4), intent(in)    :: lblock
      integer(i4), intent(in)    :: kblock
      real(r8),    intent(inout) :: farray(lblock,nnode) !inout for farray=farray_new
      integer(i4), intent(in)    :: nnode_new
      integer(i4), intent(in)    :: lblock_new
      integer(i4), intent(in)    :: kblock_new
      real(r8),    intent(inout) :: farray_new(lblock_new,nnode_new)
      type(mpi_schedule), intent(inout) :: mpi0
!--Tmps
      integer(i4) :: ii,in,iproc,ip,iip,iip2
      integer(i4) :: kptr,ncount

!-----------------------------------------------------------------
!--Allocate buffer arrays for mpi communication
       call my_allocate(mpi0%fbuff_send,mpi0%nbuff_send,'mpi_to_gpt_r8:mpi0%fbuff_send')
       call my_allocate(mpi0%fbuff_recv,mpi0%nbuff_recv,'mpi_to_gpt_r8:mpi0%fbuff_recv')
!-----------------------------------------------------------------
!------------------------------------------
!--Fill send buffers
      do ii=1,mpi0%nbuff_recv
         mpi0%fbuff_recv(ii) = farray(kblock,mpi0%ilocal_recv(ii))
      enddo
!------------------------------------------
      mpi0%itype_recv = 10000
      mpi0%itype_send = MPI_ANY_TAG

!--Post RECVs
      do iip=1,mpi0%nproc_send
        iproc  = mpi0%iproc_send(iip)-1
        kptr   = mpi0%ipntr_send(iip)
        ncount =  mpi0%ipntr_send(iip+1) - mpi0%ipntr_send(iip)
        call MPI_IRECV(mpi0%fbuff_send(kptr),ncount,MPI_REAL8,     &
                       iproc,mpi0%itype_send,MPI_COMM_WORLD,       &
                       mpi0%msgid(iip),mpi0%ierr_send)

      enddo

!--Perform SENDS
      do iip=1,mpi0%nproc_recv
        iip2   = iip + mpi0%nproc_send
        iproc  = mpi0%iproc_recv(iip)-1
        kptr   = mpi0%ipntr_recv(iip)
        ncount =  mpi0%ipntr_recv(iip+1) - mpi0%ipntr_recv(iip)
        call MPI_ISEND(mpi0%fbuff_recv(kptr),ncount,MPI_REAL8,    &
                       iproc,mpi0%itype_recv,MPI_COMM_WORLD,      &
                       mpi0%msgid(iip2),mpi0%ierr_recv)

      enddo

!--Wait for All Messages to Arrive
      call MPI_WAITALL(mpi0%nproc_all,mpi0%msgid,mpi0%istatus,mpi0%ierr)


!--Empty Recv buffers
      do ii=1,mpi0%nbuff_send
        farray_new(kblock_new,mpi0%ilocal_send(ii)) = mpi0%fbuff_send(ii)
      enddo

      call my_mpi_barrier(1) !Use barrier after waitall

!-----------------------------------------------------------------
!--Deallocate buffer arrays for mpi communication
       call my_deallocate(mpi0%fbuff_send,mpi0%nbuff_send,'mpi_to_gpt_r8:mpi0%fbuff_send')
       call my_deallocate(mpi0%fbuff_recv,mpi0%nbuff_recv,'mpi_to_gpt_r8:mpi0%fbuff_recv')
!-----------------------------------------------------------------

end subroutine mpi_to_gpt_r82
!-------------------------------------------------------------------------------
