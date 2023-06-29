!----------------------------------------------------------------------------
module mpi_schedule_read_mod
contains

  subroutine mpi_schedule_read(mpi0,iunit)

  use my_kinddefs
  use mpi_schedule_typedef
  use mp_stuff
  use io_params
  use my_allocate_mod

      type(mpi_schedule), intent(out) :: mpi0
      integer(i4),        intent(in) :: iunit
!-Tmps
      integer(i4) :: i
      integer(i4) :: idum

      idum = 0

      if (id_proc == 0) write(iwrit,*) 'Allocating and Reading in MPI Schedule ...'

      read(iunit) mpi0%nproc_send,mpi0%nproc_recv

       call my_allocate(mpi0%iproc_send,mpi0%nproc_send,'iproc_send:mpi_schedule_read')
       call my_allocate(mpi0%ipntr_send,mpi0%nproc_send+1,'iproc_send:mpi_schedule_read')
       call my_allocate(mpi0%iproc_recv,mpi0%nproc_recv,'iproc_recv:mpi_schedule_read')
       call my_allocate(mpi0%ipntr_recv,mpi0%nproc_recv+1,'iproc_recv:mpi_schedule_read')


      do i=1,mpi0%nproc_send
      read(iunit) mpi0%iproc_send(i),mpi0%ipntr_send(i)
      enddo
      read(iunit) mpi0%ipntr_send(mpi0%nproc_send+1)

      do i=1,mpi0%nproc_recv
      read(iunit) mpi0%iproc_recv(i),mpi0%ipntr_recv(i)
      enddo
      read(iunit) mpi0%ipntr_recv(mpi0%nproc_recv+1)

      mpi0%nbuff_send = mpi0%ipntr_send(mpi0%nproc_send+1) -1
      mpi0%nbuff_recv = mpi0%ipntr_recv(mpi0%nproc_recv+1) -1
      call my_allocate(mpi0%ilocal_send,mpi0%nbuff_send,'ilocal_send:mpi_schedule_read')
      call my_allocate(mpi0%ilocal_recv,mpi0%nbuff_recv,'ilocal_recv:mpi_schedule_read')


      do i=1,mpi0%nbuff_send
      read(iunit) mpi0%ilocal_send(i)
      enddo
      read(iunit) idum    !Required to match pre_nsu3d

      do i=1,mpi0%nbuff_recv
      read(iunit) mpi0%ilocal_recv(i)
      enddo
      read(iunit) idum    !Required to match pre_nsu3d

      call my_allocate(mpi0%msgid,  mpi0%nproc_send + mpi0%nproc_recv,'msgid:mpi_schedule_read')
      call my_allocate(mpi0%istatus,mpi0%nproc_send + mpi0%nproc_recv,'istatus:mpi_schedule_read')

  end subroutine mpi_schedule_read
end module mpi_schedule_read_mod
