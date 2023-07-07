!-------------------------------------------------------------------------------
module get_filename1_mod
contains

  subroutine get_filename1(output_directory,id,file_part)
  use my_kinddefs
  implicit none

      character*200, intent(in)  :: output_directory
      integer(i4),   intent(in)  :: id
      character*200, intent(out) :: file_part
!--Tmps
      character *6 ::   cpart
      integer(i4)  ::   ilen
!----------------------------------------------------------------
!--Construct Directory Name : dir_name

      if (id .lt. 10) then
      write(cpart,'(i1)') id
      elseif (id .lt. 100) then
      write(cpart,'(i2)') id
      elseif (id .lt. 1000) then
      write(cpart,'(i3)') id
      elseif (id .lt. 10000) then
      write(cpart,'(i4)') id
      elseif (id .lt. 100000) then
      write(cpart,'(i5)') id
      else
      write(cpart,'(i6)') id
      endif


      ilen = len_trim(output_directory)                !Character length (w/o blanks)
      write(file_part  ,'(a)') output_directory(1:ilen)// '/cell.' //cpart

      file_part = trim(adjustl(file_part))     !Final trimmed name

!----------------------------------------------------------------

  end subroutine get_filename1
end module get_filename1_mod
