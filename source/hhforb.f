      subroutine hhforb(f,mmxl,mxlcmp,ng,dr,xr,meshr,rstr,ncmpx,hforb
     &                 ,ls,wk)
c-----------------------------------------------------------------------
c     This routine calculates orbital hyperfine field that exist
c     spin-orbit coupling is included. The orbital hyperfine field
c     is given by
c     Hhf_orb = 2 * mu_B * sum_(l,m,s) int R*_lms (1/r**3) R_lms dr
c     Here R_lms is the radial wave functions r*Psi.
c     Numerically this is calculated as
c     Hhf_orb = 125.168 sum_(l,m,s) int R*_lms (1/r**3) R_lms dr
c     in units of G, where <1/r**3> is givein in a.u. The energy
c     integration over the occupied states has to be taken.
c     coded by H.Akai, Tokyo, 22 Jun. 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 xr(meshr,ncmpx),dr(meshr,ncmpx),rstr(meshr,mmxl,ng,ncmpx,2)
     &         ,f(mmxl,ncmpx,ng,2),hforb(ncmpx),wk(meshr,2)
      integer mxlcmp(ncmpx)
      logical ls
c     --- bmagnt: Bohr magneton in units of erg/G
c         rbohr:  Bohr radius in units of cm
      data bmagnt/9.274015d-21/, rbohr/5.29177d-9/
      if(ls) then
c     --- b=125.168
      b=2d0*bmagnt/rbohr**3*1d-3
      do 10 i=1,ncmpx
      do 20 k=1,meshr-1
   20 wk(k,2)=dr(k,i)/xr(k,i)**3
      hforb(i)=0d0
      do 10 is=1,2
      do 10 j=1,mxlcmp(i)**2
      l=lindx(j)
      m=j-l*(l-1)-1
      c=b*dble(m)
      do 10 n=1,ng
      do 30 k=1,meshr-1
   30 wk(k,1)=rstr(k,j,n,i,is)**2*wk(k,2)
   10 hforb(i)=hforb(i)+c*sintgr(wk,meshr)*f(j,i,n,is)
      else
      call clrarr(hforb,ncmpx)
      endif
      end
