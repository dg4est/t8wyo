!-------------------------------------------------------------------------------
!module cell3d_flux_test_mod
!contains


  subroutine cell3d_flux_test(ntetra,ntetra_new,ndc4,                &
                              npyr,npyr_new,ndc5,                    &
                              nprizm,nprizm_new,ndc6,                &
                              nhex,nhex_new,ndc8,                    &
                              nface3,ndf3,                           &
                              nface4,ndf4,                           &
                              nnode,xgeom,                           &
                              mpi_xcell,IMODE)
    use my_kinddefs
    use mp_stuff
    use io_params
    use mpi_schedule_typedef
    ! use mpi_from_gpt_mod
    use mpi_to_gpt_mod
    use my_allocate_mod
    use my_deallocate_mod
    use all_reduce_sum_mod
    use all_reduce_max_mod
    use triarea_mod
    implicit none

      integer(i4), intent(in)   :: ntetra,ntetra_new
      integer(i4), intent(in)   :: ndc4(4,ntetra_new)
      integer(i4), intent(in)   :: npyr,npyr_new
      integer(i4), intent(in)   :: ndc5(5,npyr_new)
      integer(i4), intent(in)   :: nprizm,nprizm_new
      integer(i4), intent(in)   :: ndc6(6,nprizm_new)
      integer(i4), intent(in)   :: nhex,nhex_new
      integer(i4), intent(in)   :: ndc8(8,nhex_new)
      integer(i4), intent(in)   :: nface3
      integer(i4), intent(inout):: ndf3(5,nface3)
      integer(i4), intent(in)   :: nface4
      integer(i4), intent(inout):: ndf4(6,nface4)
      integer(i4), intent(in)   :: nnode
      real(r8),    intent(in)   :: xgeom(3,nnode)
      type(mpi_schedule), intent(inout) :: mpi_xcell
      integer(i4), intent(in)   :: IMODE


!--Tmps
      integer(i4) :: i,n1,n2,n3,n4,ic1,ic2,ncell,ibad,ibad_all
      real(r8), pointer :: residual(:) => null()
      real(r8)    :: u1,u2,u3,flux
      real(r8)    :: fn1,fn2,fn3
      real(r8)    :: tol,fmax,fmax_all

      tol = 1.e-6

!--Test for 0 entries in face to cell pointers
      ibad  = 0
      do i=1,nface3
        if (ndf3(4,i) == 0) ibad  = ibad  + 1
        if (ndf3(5,i) == 0) ibad  = ibad  + 1
      enddo
      do i=1,nface4
        if (ndf4(5,i) == 0) ibad  = ibad  + 1
        if (ndf4(6,i) == 0) ibad  = ibad  + 1
      enddo
      call all_reduce_sum(ibad,ibad_all)
      if (ibad_all /= 0) then
        if (id_proc == 0) then
          write(iwrit,*) 'ERROR: Found 0 entries in face to cell arrays'
        endif
        call stop_all
      endif


      if (id_proc == 0) then
        write(iwrit,*) ' '
        write(iwrit,*) 'Testing face-cell connectivity for tolerance = ',tol
      endif

      ncell = ntetra_new + npyr_new + nprizm_new + nhex_new
      u1 = 1.5
      u2 = 2.5
      u3 = 5.5

      call my_allocate(residual,ncell,'test_cell_face_connectivity:residual')

        do i=1,ncell
         residual(i) = 0.0
        enddo

        do i=1,nface3
          n1 = ndf3(1,i)
          n2 = ndf3(2,i)
          n3 = ndf3(3,i)
          ic1= ndf3(4,i)
          ic2= ndf3(5,i)
          call triarea(xgeom(1,n1),xgeom(2,n1),xgeom(3,n1), &
                       xgeom(1,n2),xgeom(2,n2),xgeom(3,n2), &
                       xgeom(1,n3),xgeom(2,n3),xgeom(3,n3), &
                       fn1,fn2,fn3)
          flux = fn1*u1 + fn2*u2 + fn3*u3
          if (ic2 > 0) then
            residual(ic1) = residual(ic1) - flux
            residual(ic2) = residual(ic2) + flux
          else
            residual(ic1) = residual(ic1) - flux
          endif
        enddo

        do i=1,nface4
          n1 = ndf4(1,i)
          n2 = ndf4(2,i)
          n3 = ndf4(3,i)
          n4 = ndf4(4,i)
          ic1= ndf4(5,i)
          ic2= ndf4(6,i)
          call triarea(xgeom(1,n1),xgeom(2,n1),xgeom(3,n1), &
                       xgeom(1,n2),xgeom(2,n2),xgeom(3,n2), &
                       xgeom(1,n3),xgeom(2,n3),xgeom(3,n3), &
                       fn1,fn2,fn3)
          flux = fn1*u1 + fn2*u2 + fn3*u3
          call triarea(xgeom(1,n1),xgeom(2,n1),xgeom(3,n1), &
                       xgeom(1,n3),xgeom(2,n3),xgeom(3,n3), &
                       xgeom(1,n4),xgeom(2,n4),xgeom(3,n4), &
                       fn1,fn2,fn3)
          flux = flux + fn1*u1 + fn2*u2 + fn3*u3
          if (ic2 > 0) then
            residual(ic1) = residual(ic1) - flux
            residual(ic2) = residual(ic2) + flux
          else
            residual(ic1) = residual(ic1) - flux
          endif
        enddo


!       if (IMODE == 0) then
!         call mpi_to_gpt(ncell,residual,mpi_tetra)
!         call mpi_to_gpt(ncell,residual,mpi_pyr)
!         call mpi_to_gpt(ncell,residual,mpi_prizm)
!         call mpi_to_gpt(ncell,residual,mpi_hex)
!       else
          call mpi_to_gpt(ncell,residual,mpi_xcell)
!       endif


        fmax = 0.0
        ibad = 0
        do i=1,ncell
!       do i=1,ntetra+npyr+nprizm+nhex
          fmax = max(fmax,abs(residual(i)))
          if(abs(residual(i)) > tol) ibad = ibad + 1
        enddo

        call all_reduce_sum(ibad,ibad_all)
        call all_reduce_max(fmax,fmax_all)

        if (ibad_all == 0) then
          if (id_proc == 0) then
            write(iwrit,*) '***************************************'
            write(iwrit,*) 'Passed Face-Cell Connectivity Test'
            write(iwrit,*) 'Max residual = ',fmax_all
            write(iwrit,*) '***************************************'
          endif
        else
          if (id_proc == 0) then
            write(iwrit,*) '***************************************'
            write(iwrit,*) 'Failed Face-Cell Connectivity Test '
            write(iwrit,*) 'Number of non zero cell residuals: ',ibad_all
            write(iwrit,*) 'Max residual = ',fmax_all
            write(iwrit,*) '***************************************'
          endif
          call stop_all
        endif


      call my_deallocate(residual,ncell,'test_cell_face_connectivity:residual')


  end subroutine cell3d_flux_test
!end module test_cell_face_connectivity_mod
