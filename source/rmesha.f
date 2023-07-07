      subroutine rmesha(x1,x2,x3,dr,xr,meshr)
c-----------------------------------------------------------------------
c     +----------------------------------------------------+
c         radial mesh of the form r(x)=c*x**a*(x+g)**b+d
c     +----------------------------------------------------+
c     coded by H.Akai, 1984, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 xr(meshr),dr(meshr)
c     --- For gga exchange-correlation, a=6 could generate too fine
c         a mesh to obtain accurate 1st and 2nd derivatives of particle
c         densit, causing convergence problems. In that case, a=4 or
c         5 might give a better performance.
      data a,b,g,e/6d0,-4d0,1d-1,-10d0/
c     data a,b,g,e/5d0,-3d0,1d-1,-10d0/
      if(mod(meshr,2) .ne. 0)
     &  call errtrp(1,'rmesha','odd  meshr is illegal')
c
   20 r1=(1d0-e)/(dble(meshr-1)-e)
      c=(x1-x2)/(r1**a*(r1+g)**b-(1d0+g)**b)
      d=x2-c*(1d0+g)**b
      do 10 k=1,meshr-1
      x=(dble(k)-e)/(dble(meshr-1)-e)
      xr(k)=c*x**a*(x+g)**b+d
   10 dr(k)=(xr(k)-d)*(a/x+b/(x+g))/(dble(meshr-1)-e)
      xr(meshr)=x3
      dr(meshr)=x3-x2
c     write(6,1100)xr(1),xr(2),xr(meshr-2),xr(meshr-1)
c     write(6,1200)dr(1),dr(2),dr(meshr-2),dr(meshr-1)
c1100 format('   x1,x2,xmeshr-2,xmeshr-1=',1p,4e14.7)
c1200 format('   d1,d2,dmeshr-2,dmeshr-1=',1p,4e14.7)
      return
      end
