!***************************************************
!                 module my_kinddefs:
!     Specifies machine independant precision of 
!     double and single precision real variables 
!      and 4 types of integer variables.  The reals
!      are mad acessible through variables r4 and r8 
!      while the integers are i1, i2 i4, i8, each of 
!      which has that many digits behind the decimal place
!**********************************************************      

module my_kinddefs

  implicit none

  public r4, r8, i1, i2, i4, i8, i48
  public r4_bytesize, r8_bytesize, i1_bytesize, i2_bytesize, i4_bytesize, i8_bytesize, i48_bytesize
  public PI, TWOPI, PIO2, half, onethird, twothird, threeights, fourthird

  private
  !---> All these a based on the default IEEE 754 format 

  !---> Single Precision Real: IEEE 754
  integer, parameter :: r4=selected_real_kind(6,37)
  integer, parameter :: r4_bytesize=4
  
  !---> Double Precision Real: IEEE 754
  integer, parameter :: r8=selected_real_kind(15,307)
  integer, parameter :: r8_bytesize=8
  
  !---> Single decimal place integer (1 byte)
  integer, parameter :: i1=selected_int_kind(2)
  integer, parameter :: i1_bytesize=1
  
  !---> 2 decimal place integer (2 bytes): 
  integer, parameter :: i2=selected_int_kind(3)
  integer, parameter :: i2_bytesize=2

  !---> 4 decimal place integer ( 4 bytes):
  integer, parameter :: i4=selected_int_kind(5)
  integer, parameter :: i4_bytesize=4

  !---> 8 decimal place integer (8 bytes):
  integer, parameter :: i8=selected_int_kind(10)
  integer, parameter :: i8_bytesize=8

  !---> Variable  decimal place integer (4 or 8 bytes):
  integer, parameter :: i48=selected_int_kind(5)   !Use this for most grids...
  integer, parameter :: i48_bytesize=4
! integer, parameter :: i48=selected_int_kind(10)  !Use this for large grids on global arrays...
! integer, parameter :: i48_bytesize=8

  !---> Additionally we put some very handy constants in here
  !     to make it a bit more like matlab and other tools
  real(r8), parameter :: half = 1._r8/2._r8
  real(r8), parameter :: onethird= 1._r8/3._r8
  real(r8), parameter :: twothird= 2._r8/3._r8
  real(r8), parameter :: threeights= 3._r8/8._r8
  real(r8), parameter :: fourthird = 4._r8/3._r8
  
  real(r8), parameter :: PI=3.141592653589793238462643383279502884197_r8
  real(r8), parameter :: PIO2=1.57079632679489661923132169163975144209858_r8
  real(r8), parameter :: TWOPI=6.283185307179586476925286766559005768394_r8


end module my_kinddefs
