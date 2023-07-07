      subroutine cartes(v,n,r,ibrav)
c-----------------------------------------------------------------------
c     Given coordinates in the framework of conventional cell of any
c     lattices, this program converts them into cartesian coordinates.
c     Coded by H. Akai, Tokyo, Aug. 30, 2020
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v(3,n),r(3,3),w(3),x(3,3)
c     --- first construct conventional cell from primitive cell.
c
c     --- I (body centered) cases
      if(ibrav .eq. 2 .or. ibrav .eq. 5. .or. ibrav .eq. 8) then
      do 10 i=1,3
      x(i,1)=r(i,2)+r(i,3)
      x(i,2)=r(i,3)+r(i,1)
   10 x(i,3)=r(i,1)+r(i,2)
c
c    --- F (face centered) cases
      else if(ibrav .eq. 1 .or. ibrav .eq. 7 .or. ibrav .eq. 15) then
      do 20 i=1,3
      x(i,1)=r(i,2)+r(i,3)-r(i,1)
      x(i,2)=r(i,3)+r(i,1)-r(i,2)
   20 x(i,3)=r(i,1)+r(i,2)-r(i,3)
c
c     --- C (base centered) cases
      else if(ibrav .eq. 9 .or. ibrav .eq. 11) then
      do 30 i=1,3
      x(i,1)=r(i,1)+r(i,2)
      x(i,2)=r(i,2)-r(i,1)
   30 x(i,3)=r(i,3)
c
c     --- P (primitive) cases
      else
      do 40 j=1,3
      do 40 i=1,3
   40 x(i,j)=r(i,j)
      endif
c
c     --- now cartesian coordinates are calculated
      do 50 j=1,n
      do 60 i=1,3
      w(i)=0d0
      do 60 k=1,3
   60 w(i)=w(i)+x(i,k)*v(k,j)
      do 50 i=1,3
   50 v(i,j)=w(i)+1d-10
      end
