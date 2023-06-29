module params
  use my_kinddefs

  integer(i4) :: IERROR
! integer(i4) :: IVECTOR
! integer(i4) :: LINES_PER_GROUP

  integer(i4) :: npart
! integer(i4) :: IVGRID
! integer(i4) :: IDIST_METHOD_DEFAULT,IDIST_METHOD,IDIST_MODE
! integer(i4) :: IMAKEPOIN1
! integer(i4) :: NLEVEL_AMG
! integer(i4) :: NLAYER_REMOVE
! integer(i4) :: KLINE_COARSE
! integer(i4) :: KBFACE_PERIODIC
! integer(i4) :: KBFACE_CURVED
! integer(i4) :: KAMG_ORDER
! integer(i4) :: MAX_LINE_SIZE
! integer(i4) :: MAX_LINE_DIMENSION

  character*200 case_title
  character*200 mesh_file
  character*200 bcs_file
  character*200 output_directory

  real(r8)   :: fmach
  real(r8)   :: yangle
  real(r8)   :: zangle
  real(r8)   :: re_number

end module params


