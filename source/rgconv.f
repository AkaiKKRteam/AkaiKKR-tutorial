      subroutine rgconv(a,b,vc)
c-----------------------------------------------------------------------
c     Given three vectors a's, this program calcutates their reciprocal
c     vectors b's and the volume spanned by a's. The 1st index of a and
c     b indicates the component, i.e. x, y, and z, and the 2nd index
c     specifies one of three vectors, a1, a2, a3, etc.
c     Coded by H. Akai, Tokyo, Apr.11, 2018
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(3,3),b(3,3)
c     --- real cell volume ---
      vc=(a(2,1)*a(3,2)-a(3,1)*a(2,2))*a(1,3)
     &  +(a(3,1)*a(1,2)-a(1,1)*a(3,2))*a(2,3)
     &  +(a(1,1)*a(2,2)-a(2,1)*a(1,2))*a(3,3)
c     --- primitive reciprocal lattice vectors ---
      b(1,1)=(a(2,2)*a(3,3)-a(3,2)*a(2,3))/vc
      b(2,1)=(a(3,2)*a(1,3)-a(1,2)*a(3,3))/vc
      b(3,1)=(a(1,2)*a(2,3)-a(2,2)*a(1,3))/vc
      b(1,2)=(a(2,3)*a(3,1)-a(3,3)*a(2,1))/vc
      b(2,2)=(a(3,3)*a(1,1)-a(1,3)*a(3,1))/vc
      b(3,2)=(a(1,3)*a(2,1)-a(2,3)*a(1,1))/vc
      b(1,3)=(a(2,1)*a(3,2)-a(3,1)*a(2,2))/vc
      b(2,3)=(a(3,1)*a(1,2)-a(1,1)*a(3,2))/vc
      b(3,3)=(a(1,1)*a(2,2)-a(2,1)*a(1,2))/vc
      vc=abs(vc)
      end
