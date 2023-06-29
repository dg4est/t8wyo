!-----------------------------------------------------------------------------
subroutine start_time
  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"

  call get_time_mp(time1,time1_mp)

end subroutine start_time
!-----------------------------------------------------------------------------
subroutine report_time(ctag)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"

  character(*),  intent(in)  :: ctag
  real(r8) :: time_diff
  real(r8) :: time_diff_mp

  number_timings = number_timings + 1


  call get_time_mp(time2,time2_mp)
  time_diff    = time2 - time1
  time_diff_mp = time2_mp - time1_mp

  time0(number_timings)    = time_diff
  time0_mp(number_timings) = time_diff_mp
  ctime(number_timings)    = ctag
 
  if (id_proc == 0) then
     write(iwrit,601)
     write(iwrit,602) ctag, time_diff_mp
     write(iwrit,603)
  endif

601  format(/'--------------------------------------------------------')
602  format(a,e12.4, ' secs')
603  format('--------------------------------------------------------'/)
end subroutine report_time
!-----------------------------------------------------------------------------
subroutine report_total_time

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"

  integer(i4) :: ntimes
  real(r8) :: time_diff
  real(r8) :: time_diff_mp
  real(r8) :: fratio
  integer(i4) :: i


  call get_time_mp(time2,time2_mp)
  time_diff    = time2 - time_init
  time_diff_mp = time2_mp - time_init_mp

  if (id_proc == 0) then
     write(iwrit,601)
     write(iwrit,602) time_diff_mp
     do i=1,number_timings
       fratio = 0.0
       if (time_diff_mp > 1.e-12) then
         fratio = time0_mp(i) / time_diff_mp
       end if
       write(iwrit,612) ctime(i),time0_mp(i),fratio
     enddo
     write(iwrit,603)
  endif

601  format(/'------------------------------------------------------------------------')
602  format('TOTAL WALL CLOCK TIME:                  ',e12.4, ' secs')
612  format(a40,e12.4, ' secs;  Ratio: 'f5.2)
603  format('------------------------------------------------------------------------'/)
end subroutine report_total_time

!-----------------------------------------------------------------------------
subroutine get_time_mp(time0,time0_mp)

!--Timing Routine for Various Platforms and Multi-Processing 

  use my_kinddefs
  use io_params
  implicit none
#include "mympif.h"

  real(r8) time0,time0_mp

!--Tmps
#if defined (INTEL)
  real(r4) etime
#endif
  real(r4) tarray(2)
  real(r8) ibmtime,rtc
!
!--Timing on Single (local) Processor
!
#if defined (SV1) || (T3E)
      call second(time0)
!cdm*  call secondr(time0)
#elif defined IBM
      ibmtime = rtc()
      time0   = ibmtime
#elif defined INTEL
      time0  = etime(tarray)
      time0  = tarray(1)
#else 
      call etime(tarray)
      time0  = tarray(1)
#endif
!
!--MP Timing
!
#ifdef MPI_ON
      time0_mp  = mpi_wtime()
#else
      time0_mp  = time0
#endif

end subroutine get_time_mp
!-----------------------------------------------------------------------------

