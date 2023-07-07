      subroutine tltcry(r,g,atmicp,natm,angl)
c-----------------------------------------------------------------------
c     Given a rotation angles, this program rotate r, g, and
c     atmicp by an Euler angle (angl3, angl2, angl1).
c     Coded by H. Akai, Tokyo, Jan. 22, 2019.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 r(3,3),g(3,3),atmicp(3,natm),angl(3),rtilt(3)
      pi=4d0*atan(1d0)
      rad=pi/180d0
      angl1=rad*angl(1)
      angl2=rad*angl(2)
      angl3=rad*angl(3)
      do 10 i=1,3
      call vrotat(rtilt,r(1,i),angl3,angl2,angl1)
      do 10 j=1,3
   10 r(j,i)=rtilt(j)
      do 30 i=1,natm
      call vrotat(rtilt,atmicp(1,i),angl3,angl2,angl1)
      do 30 j=1,3
   30 atmicp(j,i)=rtilt(j)
      do 50 i=1,3
      call vrotat(rtilt,g(1,i),angl3,angl2,angl1)
      do 50 j=1,3
   50 g(j,i)=rtilt(j)
      end
