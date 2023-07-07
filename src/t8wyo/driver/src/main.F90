program main
    use my_kinddefs
    use mpi
    implicit none

    integer(i4) :: rank,nrank,comm=MPI_COMM_WORLD
    integer(i4) :: mpierr

    !comm = MPI_COMM_WORLD

    call MPI_Init(mpierr)
    call MPI_Comm_size(comm,nrank,mpierr)
    call MPI_Comm_rank(comm,rank,mpierr)

    print*,"HELLO DRIVER: rank ",rank,"of ",nrank
    call t8wyo_interface_fortran_init(comm)

    call MPI_Finalize(mpierr)
end program