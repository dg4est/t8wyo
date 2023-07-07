module my_typedefs

    use my_kinddefs
    use mpi_schedule_typedef

    type int_ptr_array
      integer(i4), pointer :: ptr(:) => null()
    end type int_ptr_array
    integer, parameter :: int_ptr_size=8

    type real_ptr_array
      real(r8)   , pointer :: ptr(:) => null()
    end type real_ptr_array
    integer, parameter :: ireal_ptr_size=8


    type mpi_array
      type(mpi_schedule) :: ptr
    end type mpi_array
    integer, parameter :: mpi_ptr_size=136



end module my_typedefs

