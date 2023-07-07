module global_arrays_dynamic
  use my_kinddefs

  integer(i48)         :: ntetra_global,npyr_global,nprizm_global,nhex_global
  integer(i48)         :: nnode_global
  integer(i48)         :: nbnode_global !(not used in mcell file)
  integer(i48)         :: nbface3_global,nbface4_global
! integer(i4), pointer :: ipart(:)     => null()
! integer(i4), pointer :: ireorder(:)  => null()
! integer(i48),pointer :: ivtxdist(:)  => null()

end module
