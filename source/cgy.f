      function cgy(l1,l2,l3,m1,m2,m3)
c-----------------------------------------------------------------------
c     Given l1,l2,l3,m1,m2,m3, this function returns the Gaunt number
c     (coefficients) int Y*_L1 Y_L2 Y_L3. Here Y_L's are the spherical
c     harmonics. The expression
c     cgy(l1,l2,l3;m1,m2,m3)=sqrt((2*l3+1)*(2*l2+1)/(2*l1+1)/4pi)
c         *c(l3,l2,l1;m3,m2,m1)*c(l3,l2,l1;0,0,0),
c     where c is the Clebsch-Gordan coefficients, is used. The
c     definition is slightly different from that given by Rose
c     (Elementary Thoury of Angular Momentum, p62) where the same
c     quantities are defined as int Y*_L3 Y_L2 Y_L1
c     
c     coded by H. Akai, 7 Feb. 2017, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      pi=4d0*atan(1d0)
      cgy=0d0
      if(l1.lt.0.or.l2.lt.0.or.l3.lt.0) return
      if(l1.lt.abs(m1).or.l2.lt.abs(m2).or.l3.lt.abs(m3)) return
      if(l3.lt.abs(l1-l2).or.l3.gt.l1+l2) return
      if(m2+m3 .ne. m1) return
      cgy=sqrt(dble((2*l3+1)*(2*l2+1))/dble(2*l1+1)/4d0/pi)
     &    *cgc(l3,l2,l1,m3,m2,m1)*cgc(l3,l2,l1,0,0,0)
      end
