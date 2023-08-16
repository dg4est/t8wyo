!----------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine cell3d_read_mesh_interface(ntetra_app,npyr_app,nprizm_app,nhex_app,             &
                                      ntetra_ng_app,npyr_ng_app,nprizm_ng_app,nhex_ng_app, &
                                      nnode_app,ngpt_app,                                  &
                                      nbface3_app,nbface4_app,nface3_app,nface4_app,       &
                                      ndc4_app,ndc5_app,ndc6_app,ndc8_app,                 &
                                      nbf3_app,nbf4_app,                                   &
                                      ifpat3_app,ifpat4_app,                               &
                                      xgeom_app,                                           &
                                      ndf3_app,ndf4_app) bind(C,name="cell3d_read_mesh_interface_")
    use iso_c_binding
    use my_kinddefs
    use mp_stuff
    use my_mpi_barrier_mod
    use local_arrays_dynamic
    implicit none

    integer(i4), intent(out) :: ntetra_app,npyr_app,nprizm_app,nhex_app
    integer(i4), intent(out) :: ntetra_ng_app,npyr_ng_app,nprizm_ng_app,nhex_ng_app
    integer(i4), intent(out) :: nnode_app,ngpt_app
    integer(i4), intent(out) :: nbface3_app,nbface4_app,nface3_app,nface4_app
    type(c_ptr),intent(out)  :: ndc4_app
    type(c_ptr),intent(out)  :: ndc5_app
    type(c_ptr),intent(out)  :: ndc6_app
    type(c_ptr),intent(out)  :: ndc8_app
    type(c_ptr),intent(out)  :: nbf3_app
    type(c_ptr),intent(out)  :: nbf4_app
    type(c_ptr),intent(out)  :: ifpat3_app
    type(c_ptr),intent(out)  :: ifpat4_app
    type(c_ptr),intent(out)  :: xgeom_app
    type(c_ptr),intent(out)  :: ndf3_app
    type(c_ptr),intent(out)  :: ndf4_app

    !--Read mesh
    if(id_proc == 0)then
        call cell3d_read_mesh()

        !--Set application mesh variables
        ntetra_app = ntetra
        npyr_app = npyr
        nprizm_app = nprizm
        nhex_app = nhex

        ntetra_ng_app = ntetra_ng
        npyr_ng_app = npyr_ng
        nprizm_ng_app = nprizm_ng
        nhex_ng_app = nhex_ng

        nnode_app = nnode
        ngpt_app = ngpt

        nbface3_app = nbface3
        nbface4_app = nbface4
        nface3_app = nface3
        nface4_app = nface4

        !--Set C pointers to Fortran pointers
        ndc4_app = c_loc(ndc4)  ! tet connectivity (4,ntet+ntet_ng)
        ndc5_app = c_loc(ndc5)  ! pyr connectivity (5,npyr+npyr_ng)
        ndc6_app = c_loc(ndc6)  ! prz connectivity (6,nprizm+nprism_ng)
        ndc8_app = c_loc(ndc8)  ! hex connectivity (8,nhex+nhex_ng)
        nbf3_app = c_loc(nbf3)  ! tri  boundary face connectivity (3,nbface3): nodes 1,2,3
        nbf4_app = c_loc(nbf4)  ! quad boundary face connectivity (4,nbface4): nodes 1,2,3,4
        ndf3_app = c_loc(ndf3)  ! tri  face connectivity (5,nface3): nodes 1,2,3,   e1,e2 (if boundary e2=negative nbf3 id)
        ndf4_app = c_loc(ndf4)  ! quad face connectivity (6,nface4): nodes 1,2,3,4, e1,e2 (if boundary e2=negative nbf4 id)
        xgeom_app = c_loc(xgeom)! node coordinates (3,nnode+ngpt)
        ifpat3_app = c_loc(ifpat3) ! tri  face patch (nbface3)
        ifpat4_app = c_loc(ifpat4) ! quad face patch (nbface4)
    endif
    call my_mpi_barrier(0)

end subroutine

subroutine cell3d_deallocate_mesh() bind(C,name="cell3d_deallocate_mesh_")
    use iso_c_binding
    use my_kinddefs
    use io_params
    use mp_stuff
    use local_arrays_dynamic
    use my_deallocate_mod
    implicit none

    if(id_proc == 0) write(iwrit,*) 'Deallocating mesh arrays...'

    if(id_proc == 0)then
        call my_deallocate(ndc4,4*(ntetra+ntetra_ng),'cell3d_allocate_mesh_data:ndc4')
        call my_deallocate(ndc5,5*(npyr+npyr_ng) ,   'cell3d_allocate_mesh_data:ndc5')
        call my_deallocate(ndc6,6*(nprizm+nprizm_ng),'cell3d_allocate_mesh_data:ndc6')
        call my_deallocate(ndc8,8*(nhex+nhex_ng) ,   'cell3d_allocate_mesh_data:ndc8')

        call my_deallocate(nbf3,3*nbface3  ,   'cell3d_allocate_mesh_data:nbf3')
        call my_deallocate(nbf4,4*nbface4  ,   'cell3d_allocate_mesh_data:nbf4')
        call my_deallocate(ifpat3,nbface3  ,   'cell3d_allocate_mesh_data:ifpat3')
        call my_deallocate(ifpat4,nbface4  ,   'cell3d_allocate_mesh_data:ifpat4')

        call my_deallocate(xgeom,3*(nnode+ngpt),'cell3d_allocate_mesh_data:xgeom')

        call my_deallocate(ndf3,5*nface3,   'cell3d_allocate_mesh_data:ndf3')
        call my_deallocate(ndf4,6*nface4,   'cell3d_allocate_mesh_data:ndf4')
    endif
end subroutine

subroutine cell3d_read_mesh()
    use my_kinddefs
    use mp_stuff
    use my_mpi_barrier_mod
    use get_filename1_mod
    use io_params
    use params,only:mesh_file
    use global_arrays_dynamic
    use local_arrays_dynamic
    use mpi_schedules
    implicit none
    character*200 file_part

    !--Construct partition file name to be read in
    call get_filename1(mesh_file,id_proc,file_part)

    !--Open partition mesh file
    open(unit=io_mesh,file=file_part,form='unformatted',access='STREAM')

    call cell3d_read_mesh_size()
    call cell3d_allocate_mesh_data()
    call cell3d_read_mesh_data(ntetra,npyr,nprizm,nhex,             &
                               ntetra_ng,npyr_ng,nprizm_ng,nhex_ng, &
                               nnode,ngpt,                          &
                               nbface3,nbface4,nface3,nface4,       &
                               ndc4,ndc5,ndc6,ndc8,                 &
                               nbf3,nbf4,                           &
                               ifpat3,ifpat4,                       &
                               xgeom,                               &
                               mpi_xcell,mpi_gpt,                   &
                               ndf3,ndf4)
    close(io_mesh)

    if(id_proc == 0) then
        write(6,*) 'DONE *****************'
    endif

end subroutine cell3d_read_mesh

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine cell3d_read_mesh_size
  use io_params
  use mp_stuff
  use global_arrays_dynamic
  use local_arrays_dynamic
  implicit none

  if(id_proc == 0) write(iwrit,*) 'Reading in mesh size params...'

!--READ in global mesh sizes
      read(io_mesh) ntetra_global,npyr_global,nprizm_global,nhex_global
      read(io_mesh) nnode_global,nbface3_global,nbface4_global,npatch,ncomp,nbod

!--READ in local partition sizes
      read(io_mesh) ntetra,npyr,nprizm,nhex,nbface3,nbface4,nnode
      read(io_mesh) ntetra_ng,npyr_ng,nprizm_ng,nhex_ng,ngpt
      read(io_mesh) nface3,nface4

      if(id_proc == 0) write(iwrit,601) ntetra_global,npyr_global,nprizm_global,nhex_global,nnode_global, &
                                        nbface3_global,nbface4_global,npatch,ncomp,nbod
601   format(/'----------------------------------------------------------------------------',/&
              '                       MESH GLOBAL STATISTICS',/                               &
              '----------------------------------------------------------------------------',/&
              '    NTETRA      NPYR    NPRIZM      NHEX     NNODE',/,5i10,/, &
              '   NBFACE3   NBFACE4    NPATCH     NCOMP      NBOD',/,5i10,/,&
              '----------------------------------------------------------------------------',/)

end subroutine cell3d_read_mesh_size
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine cell3d_allocate_mesh_data
  use io_params
  use mp_stuff
  use local_arrays_dynamic
  use my_allocate_mod
  implicit none

  if(id_proc == 0) write(iwrit,*) 'Allocating mesh arrays...'

        call my_allocate(ndc4,4*(ntetra+ntetra_ng),'cell3d_allocate_mesh_data:ndc4')
        call my_allocate(ndc5,5*(npyr+npyr_ng) ,   'cell3d_allocate_mesh_data:ndc5')
        call my_allocate(ndc6,6*(nprizm+nprizm_ng),'cell3d_allocate_mesh_data:ndc6')
        call my_allocate(ndc8,8*(nhex+nhex_ng) ,   'cell3d_allocate_mesh_data:ndc8')

        call my_allocate(nbf3,3*nbface3  ,   'cell3d_allocate_mesh_data:nbf3')
        call my_allocate(nbf4,4*nbface4  ,   'cell3d_allocate_mesh_data:nbf4')
        call my_allocate(ifpat3,nbface3  ,   'cell3d_allocate_mesh_data:ifpat3')
        call my_allocate(ifpat4,nbface4  ,   'cell3d_allocate_mesh_data:ifpat4')

        call my_allocate(xgeom,3*(nnode+ngpt),'cell3d_allocate_mesh_data:xgeom')

        call my_allocate(ndf3,5*nface3,   'cell3d_allocate_mesh_data:ndf3')
        call my_allocate(ndf4,6*nface4,   'cell3d_allocate_mesh_data:ndf4')

end subroutine cell3d_allocate_mesh_data
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine cell3d_read_mesh_data(ntetra,npyr,nprizm,nhex,             &
                                 ntetra_ng,npyr_ng,nprizm_ng,nhex_ng, &
                                 nnode,ngpt,                          &
                                 nbface3,nbface4,nface3,nface4,       &
                                 ndc4,ndc5,ndc6,ndc8,                 &
                                 nbf3,nbf4,                           &
                                 ifpat3,ifpat4,                       &
                                 xgeom,                               &
                                 mpi_xcell,mpi_gpt,                   &
                                 ndf3,ndf4)
  use my_kinddefs
  use io_params
  use mp_stuff
  use mpi_schedule_typedef
  use mpi_schedule_read_mod

! use local_arrays_dynamic
  implicit none
  integer(i4), intent(in) :: ntetra,npyr,nprizm,nhex
  integer(i4), intent(in) :: ntetra_ng,npyr_ng,nprizm_ng,nhex_ng
  integer(i4), intent(in) :: nnode,ngpt
  integer(i4), intent(in) :: nbface3,nbface4,nface3,nface4
  integer(i4),intent(out) :: ndc4(4,ntetra+ntetra_ng)
  integer(i4),intent(out) :: ndc5(5,npyr+npyr_ng)
  integer(i4),intent(out) :: ndc6(6,nprizm+nprizm_ng)
  integer(i4),intent(out) :: ndc8(8,nhex+nhex_ng)
  integer(i4),intent(out) :: nbf3(3,nbface3)
  integer(i4),intent(out) :: nbf4(4,nbface4)
  integer(i4),intent(out) :: ifpat3(nbface3)
  integer(i4),intent(out) :: ifpat4(nbface4)
  real(r8)   ,intent(out) :: xgeom(3,nnode+ngpt)
  integer(i4),intent(out) :: ndf3(5,nface3)
  integer(i4),intent(out) :: ndf4(6,nface4)
  type(mpi_schedule), intent(out) :: mpi_xcell
  type(mpi_schedule), intent(out) :: mpi_gpt

  integer(i4)  :: i,k

  if(id_proc == 0) write(iwrit,*) 'Reading in mesh arrays...'

!--Read in real cells
      read(io_mesh) ((ndc4(k,i),k=1,4),i=1,ntetra)
      read(io_mesh) ((ndc5(k,i),k=1,5),i=1,npyr)
      read(io_mesh) ((ndc6(k,i),k=1,6),i=1,nprizm)
      read(io_mesh) ((ndc8(k,i),k=1,8),i=1,nhex)
!--Read in ghost cells
      read(io_mesh) ((ndc4(k,i),k=1,4),i=ntetra+1,ntetra+ntetra_ng)
      read(io_mesh) ((ndc5(k,i),k=1,5),i=npyr+1,npyr+npyr_ng)
      read(io_mesh) ((ndc6(k,i),k=1,6),i=nprizm+1,nprizm+nprizm_ng)
      read(io_mesh) ((ndc8(k,i),k=1,8),i=nhex+1,nhex+nhex_ng)

!--Read in boundary faces
      read(io_mesh) ((nbf3(k,i),k=1,3),i=1,nbface3)
      read(io_mesh) ((nbf4(k,i),k=1,4),i=1,nbface4)
      read(io_mesh) (ifpat3(i),i=1,nbface3)
      read(io_mesh) (ifpat4(i),i=1,nbface4)

!--Read in vertex coordinates
      read(io_mesh) ((xgeom(k,i),k=1,3),i=1,nnode+ngpt)

!--Read in mpi cell to cell schedule
      call mpi_schedule_read(mpi_xcell,io_mesh)
!--Read in mpi node to node schedule
      call mpi_schedule_read(mpi_gpt,io_mesh)

!--Read in faces
      read(io_mesh) ((ndf3(k,i),k=1,5),i=1,nface3)
      read(io_mesh) ((ndf4(k,i),k=1,6),i=1,nface4)

end subroutine cell3d_read_mesh_data
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
