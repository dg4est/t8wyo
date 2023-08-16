subroutine cell3d_test_connectivity

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
  integer(i4) :: ntetra_new
  integer(i4) :: npyr_new
  integer(i4) :: nprizm_new
  integer(i4) :: nhex_new
  integer(i4) :: IMODE    



        ntetra_new = ntetra + ntetra_ng
        npyr_new   = npyr   + npyr_ng
        nprizm_new = nprizm + nprizm_ng
        nhex_new   = nhex   + nhex_ng
        nnode      = nnode  + ngpt

        call cell3d_flux_test(ntetra,ntetra_new,ndc4,                &
                              npyr,npyr_new,ndc5,                    &
                              nprizm,nprizm_new,ndc6,                &
                              nhex,nhex_new,ndc8,                    &
                              nface3,ndf3,                           &
                              nface4,ndf4,                           &
                              nnode,xgeom,                           &
                              mpi_xcell,IMODE)


end subroutine cell3d_test_connectivity


