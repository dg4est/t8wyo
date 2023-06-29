module mpi_schedule_typedef

    use my_kinddefs
    type mpi_schedule

      integer(i4)          :: nproc_send
      integer(i4)          :: nproc_recv
      integer(i4)          :: nbuff_send
      integer(i4)          :: nbuff_recv

      integer(i4), pointer :: iproc_send(:)  => null()
      integer(i4), pointer :: iproc_recv(:)  => null()
      integer(i4), pointer :: ipntr_send(:)  => null()
      integer(i4), pointer :: ipntr_recv(:)  => null()

      integer(i4), pointer :: ibuff_send(:)  => null()
      integer(i4), pointer :: ibuff_recv(:)  => null()
      real(r8),    pointer :: fbuff_send(:)  => null()
      real(r8),    pointer :: fbuff_recv(:)  => null()

      integer(i4), pointer :: ilocal_send(:) => null()
      integer(i4), pointer :: ilocal_recv(:) => null()

      integer(i4), pointer ::  msgid(:)      => null()
      integer(i4), pointer ::  istatus(:)    => null()

      integer(i4)          :: ierr_send,ierr_recv,ierr
      integer(i4)          :: itype_send,itype_recv
      integer(i4)          :: nproc_all

    end type mpi_schedule

end module mpi_schedule_typedef

