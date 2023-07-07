!----------------------------------------------------------------------------
subroutine char_adjust(filen)

!--Subroutine to adjust character string filen
!--Clip all characters after first blank encounter:

  use my_kinddefs
  implicit none

  character *200 filen

!--Tmps
  character *205 filen0
  INTEGER i,nchar

      filen = trim(adjustl(filen))

      do i=1,200
      filen0(i:i) = filen(i:i)
      filen(i:i)  = ' '
      enddo

      nchar = index(filen0,' ')-1
      read(filen0(1:nchar),'(a200)') filen

end subroutine char_adjust
