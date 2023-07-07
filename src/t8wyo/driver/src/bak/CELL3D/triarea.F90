!-----------------------------------------------------------------------------
module triarea_mod
contains

      subroutine triarea(x1,y1,z1,x2,y2,z2,x3,y3,z3,area1,area2,area3)

  use my_kinddefs
  implicit none

      real(r8), intent(in)  ::  x1,y1,z1
      real(r8), intent(in)  ::  x2,y2,z2
      real(r8), intent(in)  ::  x3,y3,z3
      real(r8), intent(out) ::  area1,area2,area3
!--Tmps
      real(r8) ::  dx1,dx2,dx3
      real(r8) ::  dy1,dy2,dy3
      real(r8) ::  dz1,dz2,dz3
      real(r8) ::  xx1,xx2,xx3
      real(r8) ::  yy1,yy2,yy3
      real(r8) ::  zz1,zz2,zz3
!
!--X projection
      dz1      = z2 - z1
      dz2      = z3 - z2
      dz3      = z1 - z3
      yy1      = y2 + y1
      yy2      = y3 + y2
      yy3      = y1 + y3
      area1    = 0.5 * (yy1*dz1 + yy2*dz2 + yy3*dz3)
!--Y projection
      dx1      = x2 - x1
      dx2      = x3 - x2
      dx3      = x1 - x3
      zz1      = z2 + z1
      zz2      = z3 + z2
      zz3      = z1 + z3
      area2    = 0.5 * (zz1*dx1 + zz2*dx2 + zz3*dx3)
!--Z projection
      dy1      = y2 - y1
      dy2      = y3 - y2
      dy3      = y1 - y3
      xx1      = x2 + x1
      xx2      = x3 + x2
      xx3      = x1 + x3
      area3    = 0.5 * (xx1*dy1 + xx2*dy2 + xx3*dy3)

end subroutine triarea
end module triarea_mod
