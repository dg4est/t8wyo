
module  io_params

  use my_kinddefs
  integer(i4) :: iread,iwrit,iterm1,iterm2,iunbuf
  integer(i4) :: io_mesh 
  integer(i4) :: io_local
  integer(i4) :: io_amg
  integer(i4) :: io_poin1
  integer(i4) :: io_distf
  integer(i4) :: io_bcs
  integer(i4) :: io_comp
  integer(i4) :: io_xbface
  integer(i4) :: io_mcell
  integer(i4) :: io_mem,io_mem_proc   !For memory allocation log...
  integer(i4) :: IPARALLEL_READ
  integer(i4) :: IPARALLEL_WRITE
  integer(i4) :: IPARALLEL_AMG_WRITE
  integer(i4) :: IO_SYSTEMCALL,IO_SLEEP,IO_UNITS(100)
  integer(i4) :: IO_VERSION
  integer(i4) :: IVERBOSE
  integer(i4) :: POINTER_ASSOCIATION_CHECK = 1
  integer(i4) :: POINTER_DEASSOCIATION_CHECK = 1

end module 
