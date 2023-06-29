!----------------------------------------------------------------------
subroutine cell3d_initialize
  use my_kinddefs
  implicit none

  character(200)  :: fileinput
  fileinput = "input.1"

  call program_mpi_init
  !call readarg(fileinput)   !Omit this if running from Python Interface
  print*,"FILE READ:",fileinput
  call cell3d_init(fileinput)

end subroutine cell3d_initialize