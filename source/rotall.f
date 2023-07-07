      subroutine rotall(r,rpt,nrpt,atmicp,natm,cornr,nc,g,gpt,ngpt,angl)
c-----------------------------------------------------------------------
c     Given a rotation Euler angle, this program rotate r, rpt, atmicp
c     cornr, and gpt by the angle.
c     Coded by H. Akai, Tokyo, 2016
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 r(3,3),rpt(3,nrpt),atmicp(3,natm),cornr(3,nc),g(3,3)
     &      ,gpt(3,ngpt),angl(3),rtilt(3)
      pi=4d0*atan(1d0)
      rad=pi/180d0
      angl1=rad*angl(1)
      angl2=rad*angl(2)
      angl3=rad*angl(3)
      do 10 i=1,3
      call vrotat(rtilt,r(1,i),angl3,angl2,angl1)
      do 10 j=1,3
   10 r(j,i)=rtilt(j)
      do 20 i=1,nrpt
      call vrotat(rtilt,rpt(1,i),angl3,angl2,angl1)
      do 20 j=1,3
   20 rpt(j,i)=rtilt(j)
      do 30 i=1,natm
      call vrotat(rtilt,atmicp(1,i),angl3,angl2,angl1)
      do 30 j=1,3
   30 atmicp(j,i)=rtilt(j)
      do 40 i=1,nc
      call vrotat(rtilt,cornr(1,i),angl3,angl2,angl1)
      do 40 j=1,3
   40 cornr(j,i)=rtilt(j)
      do 50 i=1,3
      call vrotat(rtilt,g(1,i),angl3,angl2,angl1)
      do 50 j=1,3
   50 g(j,i)=rtilt(j)
      do 60 i=1,ngpt
      call vrotat(rtilt,gpt(1,i),angl3,angl2,angl1)
      do 60 j=1,3
   60 gpt(j,i)=rtilt(j)
      end
