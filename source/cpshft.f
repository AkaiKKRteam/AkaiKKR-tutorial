      subroutine cpshft(e,tchc,tchs,cm,fcs,gs,tc,mxl,mxlcmp,ng,t,pf
     &          ,str,sra)
c-----------------------------------------------------------------------
c     In this version, pf does not include the phase factor i^(l1-l2)
c     any more. It is now attached to the scattering path operator
c     calculated  by kkrsed.
c     Modified by H. Akai, 30 July 2005, Osaka
c
c     ------------------------------------
c     --- KKR-CPA + spin-orbit version ---
c     ------------------------------------
c     Given Tchebycheff expansion coeficients tchc,tchs and the
c     shift constants fcs, this program returns t-matrix (or its
c     inverse) and the phase factor fanction, t and pf.
c     To facilitate the calculation of the phase of the Green's
c     function this program also returns function 'str'.
c     coded by H.Akai, 1983, Juelich
c     Only the modification is that now 'str' is defined differently.
c     KKR CPA implemented by H.Akai, 21 Sep 1996, Osaka
c     adapted to spin-orbit version, 20 April, 1997, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c     f resonance is taken into account
c     by H. Akai, Aug. 2014, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 tc(ng),t(mxl**2),pf(mxl**2,mxl**2),str(mxl**2)
     &          ,p(36),gs(mxl**2),e,ek,ekl,ea,erel,eef,ex,cx,sx,fx,gam
      real*8 tchc(ng,mxl**2),tchs(ng,mxl**2),fcs(4,mxl**2),cm(ng,mxl**2)
      logical sra
      data c/274.0720442d0/
      mmxl=mxlcmp**2
      erel=e
      if(sra) erel=e+(e/c)**2
      ek=sqrt(erel)
      do 10 l=1,mmxl
      ll=sqrt(dble(l)-5d-1)+1
      eef=e-fcs(3,l)
      fx=(0d0,0d0)
      cx=(0d0,0d0)
      sx=(0d0,0d0)
      do 20 n=1,ng
      fx=fx+cm(n,l)*tc(n)
      cx=cx+tchc(n,l)*tc(n)
   20 sx=sx+tchs(n,l)*tc(n)
      if(abs(fcs(4,l)) .lt. 1d-6) then
      cx=cx*eef+fcs(1,l)
      sx=sx*eef+fcs(2,l)
      else
      gam=fcs(4,l)/(eef**2+fcs(4,l)**2)
      cx=cx*eef*gam+fcs(1,l)
      sx=sx*eef*gam+fcs(2,l)
      endif
      ekl=ek**(ll-1)
      ea=erel**(ll-1)
      ex=cx+(0d0,1d0)*sx*ea*ek
      p(l)=ekl/ex 
c     if(l .eq. 10) write(*,'(5f15.7)')dble(e),dble(cx),dble(sx)
c     ---str(l) is the quantity needed to calculate
c        phase of the KKR Green's function.
      str(l)=log(sx*ea)
c     ---gs(l) is the single potential Green's function.
c     gs(l)=(fx-(0d0,1d0)*ek*ekl*p(l))/cx
      gs(l)=(fx-(0d0,1d0)*ea*ek/ex)/cx
c     gs=(0d0,0d0)
c     t(l)=1d0/(cx/(sx*ea)-(0d0,1d0)*ek)
c     --- t(l) is the negative inverse of t-matrix
c     cnorm=(cx/sx)**2+erel
   10 t(l)=-ex/sx/ea
c     --- check for spin-orbit coupling
c        The phase factor i**(l2-l3) is not needed since it is
c        attached in the subroutine kkrsed. In the
c        definitioan of the Green function it should be
c        attached though omitting it is rather harmless untill
c        the very end where the phase of the wave fucntion
c        is considered.
      do 30 l1=1,mmxl
c     ll1=sqrt(dble(l1)-5d-1)+1
      do 30 l2=1,mmxl
c     ll2=sqrt(dble(l2)-5d-1)+1
   30 pf(l1,l2)=p(l1)*p(l2)
c  30 pf(l1,l2)=p(l1)*p(l2)*(0d0,1d0)**(ll1-ll2)
      return
      end
