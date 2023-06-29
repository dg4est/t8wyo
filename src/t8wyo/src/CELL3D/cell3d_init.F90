!------------------------------------------------------------------------------
subroutine cell3d_init(fileinput)

  use my_kinddefs
! use global_arrays_dynamic
  use io_params
  use params
  use mp_stuff
! use my_allocate_mod

  implicit none
#include "mympif.h"

  character*200 fileinput

!------------------------------------------------------------------------
!------------------INITIALIZE PROGRAM MODE PARAMETERS--------------------
!------------------------------------------------------------------------

  call init_io
! call pointer_nullify ! Shoudl be able to remove this when using nullify statements in modules

!-----------------------------------------------------------
!------------------------------------------------------------------------
!--------------INITIALIZE CONSTANTS--------------------------------------
!------------------------------------------------------------------------
!
!
!     call make_cellpatterns
!     call fixed_arrays_allocate
!
!     if (IDYNAMIC_MEM .eq. 1) then
!     call all_allocate0   ! Allocate fixed size arrays prior
!     endif                ! to reading in mesh size
!
!------------------------------------------------------------
!------------------READ IN PARAMETER FILE--------------------
!------------------------------------------------------------
!------------------------------------------------------------
!
      call readpm(fileinput,       &
                  case_title,      &
                  mesh_file,       &
                  bcs_file,        &
                  fmach,yangle,zangle,re_number)
!
!------------------------------------------------------------
!--Determine if AMG3d will need to be run sequentially
!     call amg_set_mode(mg_file,NLEVEL_AMG,IRUN_AMG3D)

!------------------------------------------------------------
!--Default is to use npart = mpi num_proc specified on command line
!--although non-zero npart value <= num_proc is supported as well
!
      if (npart <= 0) then
        npart = num_proc
      endif
      npart = min(npart,num_proc)

!------------------------------------------------------------
!--Allocate fixed size arrays (use num_proc)
!     call my_allocate(ivtxdist,num_proc+1,'pre_nsu3d_init:ivtxdist')
!------------------------------------------------------------
!--Output Message...
!
      call wout(case_title,      &
                mesh_file,       &
                bcs_file,        &
                fmach,yangle,zangle,re_number)

end subroutine cell3d_init
