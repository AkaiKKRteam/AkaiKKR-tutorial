      subroutine getrtm(u,alpha,beta,gamma,mxl)
c-----------------------------------------------------------------------
c     calculate Wigner's D matrix, i.e. rotation matrix for spherical
c     harmonics rotated by an Euler angle (gamma, beta, alpha),
c     following the original foumula.
c     
c     coded by H.Akai, 14 Dec. 1996, Osaka
c     modified by H. Akai Aug. 2013 based on getrtm3 coded by M. Ogura.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mmxl=7)
      complex*16 u((2*mxl-1)**2,mxl),cunit,pf1,pf2
      real*8 fact(0:2*mmxl)
      logical init
      save init,fact
      data init/.true./, cunit/(0d0,1d0)/, zero/1d-8/
      pi=4d0*atan(1d0)
      if(mxl .gt. mmxl) call errtrp(1,'getrtm','mxl too large')
c
      if(init) then
      init=.false.
      fact(0)=1d0
      do 30 i=1,2*mmxl
   30 fact(i)=fact(i-1)*dble(i)
      endif
c
      cosb=cos(5d-1*beta)
      sinb=-sin(5d-1*beta)
      do 10 l=0,mxl-1
      do 10 m2=-l,l
      pf2=exp(-cunit*gamma*dble(m2))
      do 10 m1=-l,l
      pf1=exp(-cunit*alpha*dble(m1))
c
      if(abs(beta) .lt. zero) then
      if(m1 .eq. m2) then
      dbet=1d0
      else
      dbet=0d0
      endif
      else if(abs(beta-pi) .lt. zero) then
      if(m1 .eq. -m2) then
      dbet=(-1d0)**(l+m1)
      else
      dbet=0d0
      endif
      else
      kmn=max(0,m1-m2)
      kmx=min(l+m1,l-m2)
      dbet=0d0
      xa=sqrt(fact(l+m1)*fact(l-m1)*fact(l+m2)*fact(l-m2))
      do 20 k=kmn,kmx
      xb=fact(l+m1-k)*fact(l-m2-k)*fact(k-m1+m2)*fact(k)
 20   dbet=dbet+(-1d0)**k*xa/xb*cosb**(2*l+m1-m2-2*k)*sinb**(2*k-m1+m2)
      endif
      m21=(2*l+1)*(l+m1)+l+m2+1
 10   u(m21,l+1)=pf1*dbet*pf2
      end
