module mpi_schedules
   use mpi_schedule_typedef
   type(mpi_schedule), target :: mpi_xcell
   type(mpi_schedule), target :: mpi_gpt
!  type(mpi_schedule), target :: mpi_cell
!  type(mpi_schedule), target :: mpi_edge
!  type(mpi_schedule), target :: mpi_bgpt
end module mpi_schedules
