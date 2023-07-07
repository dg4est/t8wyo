!-------------------------------------------------------------------------------
subroutine wout(case_title,      &
                mesh_file,       &
                bcs_file,        &
                fmach,yangle,zangle,re_number)

  use my_kinddefs
  use io_params
  use mp_stuff
  implicit none
#include "mympif.h"

  character*200 case_title
  character*200 mesh_file
  character*200 bcs_file
  real(r8)    :: fmach
  real(r8)    :: yangle
  real(r8)    :: zangle
  real(r8)    :: re_number
!--TMPS
  character*80 ptitle(20)
  integer(i4) :: npline
  integer(i4) :: i

!----------------------------------------------------------------------
      if (id_proc .eq. 0) then

      npline    = 14
      ptitle(1) = '           PROGRAM  CELL3D               '
      ptitle(2) = '           (VERSION 1.0.a:r-000)         '
      ptitle(3) = '              (June  2023)               '
      ptitle(4) = '                                         '
      ptitle(5) = '         CELL CENTERED TEST CODE         '
      ptitle(6) = '                  FOR                    '
      ptitle(7) = '          T8CODE (DYNAMIC AMR)           '
      ptitle(8) = '                                         '
      ptitle(9) = '                                         '
      ptitle(10)= '                                         '
      ptitle(11)= '                                         '
      ptitle(12)= 'COPYRIGHT (C) SCIENTIFIC SIMULATIONS(2023)'
      ptitle(13)= '                                         '
      ptitle(14)= '                                         '

      write(iwrit,101)
      do 100 i=1,npline
      write(iwrit,102) ptitle(i)
  100 continue
      write(iwrit,103)
  101 format(/'--------------------------------------------------------------------------------',   &
             /'********************************************************************************'/)
  102 format(18x,a50)
  103 format(/'********************************************************************************',   &
             /'--------------------------------------------------------------------------------'/)

      write (iwrit,524) trim(adjustl(case_title))

!--Write out all Parameters

      write(iwrit,501) trim(adjustl(mesh_file))
      write(iwrit,505) trim(adjustl(bcs_file))
      write(iwrit,595) fmach,yangle,zangle,re_number



      write(iwrit,104)
  104 format(/'--------------------------------------------------------------------------------'/)

    
      endif
!----------------------------------------------------------------------
      RETURN
!----------------------------------------------------------------------

  524 format(a)
  501 format('MESH CELL FILE:',/,a)
  502 format('POIN1 FILE     ',/,a)
  503 format('MESH DIST FILE:',/,a)
  504 format('AMG FILE:      ',/,a)
  513 format('***SEQUENTIAL AMG3D WILL BE RUN (amg file does not exist)***')
  505 format('BCS FILE:      ',/,a)
  506 format('COMP FILE:     ',/,a)
  509 format('NLEVELS ((MULTIGRID: AMG levels + 1):  ',/,i4)
  510 format('NPART:         ',/,i7)
  520 format('NLAYER OF PRISMS TO BE REMOVED ',/,i4)
  581 format('NUMBER CYCLES FOR EIKONAL SOLVER (optional input):',i6)
  582 format('RMS TOLERANCE FOR EIKONAL SOLVER (optional input):',e12.4)
  583 format('MAX TOLERANCE FOR EIKONAL SOLVER (optional input):',e12.4)
  595 format('      MACH    YANGLE    ZANGLE     RE_NUMBER',/,3f10.3,e15.4)

end subroutine wout
