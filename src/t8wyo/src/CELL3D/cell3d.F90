!----------------------------------------------------------------------
program cell3d

  use my_kinddefs

  call cell3d_initialize
  call cell3d_read_mesh

  call cell3d_test_connectivity

  call stop_all

end program cell3d
!----------------------------------------------------------------------

