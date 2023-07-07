      subroutine fczero(cm,ew,ez,tchc,tchs,fcs,tm,elvl
     &                 ,ng,jmx,jmxcmp,ls)
c-----------------------------------------------------------------------
c     This program searches all zeros of C and calculates the
c     contribution of the suprious poles of the pseudo Green's function.
c     The pseudo Green's function is defined as
c
c     G(E)=( -E**l/(C**2+E**(2*l+1)*S**2) * ((S/C)*E**(l+1)+i*sqrt(E))
c            + f/C ) * R(E)*R(E)
c         =(  i*E**(l+1/2)/(C+i*E**(l+1/2)*S) + f )*(1/C) * R(E)*R(E),
c
c     where l=0,1,2,..., and f(E) is a polynomial which fits 1/S(E)
c     at all E=En's satisfing C(En)=0. In the output 'cm' is the
c     Chebyshev expansion coefficients of f(E).
c     coded by H.Akai, 1983, Juelich
c     major modification by H.Akai, 1993, Osaka
c     very minor modification by H. Akai, 25 Aug. 1999, Osaka
c     A separate treatment for the f resonance states is introduced
c     by H. Akai, Aug. 2014, Tokyo
c     Modified by H. Akai, Tokyo, Jan. 2018.
c     Major revision by H. Akai, Tokyo, Oct. 12, 2020
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 tchc(ng,jmx),tchs(ng,jmx),fcs(4,jmx),cm(ng,jmx)
     &         ,tm(ng,ng),elvl(ng)
      real*8,allocatable::c(:),d(:),extr(:)
      complex*16 ce,cpolin
      complex*16,allocatable::roots(:),tc(:),y(:)
      logical ls
      data zero/1d-10/, small/1d-7/, esmear/1d-2/, elim/1.1d0/
      allocate(roots(ng),tc(ng),y(ng),d(ng+1),c(ng),extr(ng))
      nr=0
      do 10 i=1,ng
      do 10 l=1,jmx
   10 cm(i,l)=0d0
      do 20 j=1,jmxcmp
      jj=sqrt(dble(j)-5d-1)+1
      if(ls .or. j .eq. (jj-1)**2+1) then
c
c     --- polinomial expansion of c-function where
c         d's are the polinomial coefficients:
c         c=d(n)*x**(n-1)+ ... +d(1)
      call chebpc(tchc(1,j),d,ng)
c     ---working on scaled energy centered on 'ew', scaled by 'ez'.
c     --- multiply the polinomial by e-e0:
c     --- new polinomials =  (e-e0)*(Sum d(n)*e**(n-1) + ... +d(1))
      e0=(fcs(3,j)-ew)/ez
      d(ng+1)=d(ng)
      do 30 i=ng,2,-1
   30 d(i)=d(i-1)-d(i)*e0
      d(1)=-d(1)*e0
c     --- for f-states, a bit different fitting has been made
c         and additional term is needed.
c     --- First, fit the function sqrt(...) into Tchebyshev
c        polinomials, which is possible because the argument
c        of sqrt(...) never takes zero, then transform it
c        into a polinomial.
c
      if(abs(fcs(4,j)) .gt. 0d0) then
      delt=fcs(4,j)/ez
      d(1)=d(1)+fcs(1,j)*(delt+e0**2/delt)
      d(2)=d(2)-2d0*fcs(1,j)*e0/delt
      d(3)=d(3)+fcs(1,j)/delt
      else
c
c     call clrarr(extr,ng)
c     do 40 k=1,ng
c     c0=fcs(1,j)*sqrt((elvl(k)-fcs(3,j))**2+fcs(4,j)**2)/ez
c     c0=fcs(1,j)*((elvl(k)-fcs(3,j))**2+fcs(4,j)**2)/ez/fcs(4,j)
c     do 40 n=1,ng
c  40 extr(n)=extr(n)+tm(n,k)*c0
c     call chebpc(extr,c,ng)
c     do 50 i=1,ng
c  50 d(i)=d(i)+c(i)
c     else
      d(1)=d(1)+fcs(1,j)/ez
      endif
c     -- reduce the degree of the polinomial.
      big=zero
      nr0=0
      do 40 i=1,ng+1
      big=max(big,abs(d(i)))
      if(abs(d(i)/big) .gt. zero) then
      nr0=i-1
      else
      d(i)=0d0
      endif
   40 continue
c
c     --- try to find zeros of c-function
      itrf=nr0
      do 50 itry=nr0,2,-1
      nr=itry
      call zroots(d,nr,roots,ierr)
      if(ierr .eq. 0) exit
   50 itrf=itry
c     write(*,*)j,ew,ez
c     write(*,'(i3,2f12.5)')(i,ew+ez*roots(i),i=1,nr)
      if(itrf .eq. 2) call errtrp(1,'fczero','root finding fails')
      do 60 i=2,nr
      if(abs(roots(i)-roots(i-1)) .lt. small) then
c     write(*,'(1x,i3,2f15.7,3x,2f15.7)')i,roots(i),roots(i-1)
      call errtrp(1,'fczero','same zero counted twice')
      endif
   60 continue
c
c     --- Pick up only poles that lie within a circle of
c         radius elim*ez centered on ew using Fermi's distribution
c         function fd with a smearing factor esmear.
      do 70 i=1,nr
      arg=(abs(roots(i))-elim)/esmear
      arg=max(min(46d0,arg),-46d0)
      fd=1d0/(1d0+exp(arg))
c     ---now scale them back.
      roots(i)=ew+roots(i)*ez
c
c     --- y(i)=1/s(i) at the i-th zero of c-function.
      call cgntcs(roots(i),ew,ez,tc,ng)
      y(i)=(0d0,0d0)
      do 80 l=1,ng
   80 y(i)=y(i)+tchs(l,j)*tc(l)
      y(i)=y(i)*(roots(i)-fcs(3,j))
c     --- for states with a resonance
      if(abs(fcs(4,j)) .gt. 0d0)
     &     y(i)=y(i)*fcs(4,j)/((roots(i)-fcs(3,j))**2+fcs(4,j)**2)
   70 y(i)=fd/(y(i)+fcs(2,j))
c
c     --- Fit y(i) by a polinomial, and then exapnd it into
c         Tchebyshev polinomials.
      do 90 k=1,ng
      ce=dcmplx(elvl(k),0d0)
      f=dble(cpolin(ce,roots,y,nr,err))
      do 90 l=1,ng
   90 cm(l,j)=cm(l,j)+f*tm(l,k)
c     if(j .eq. 10)
c      --- check the accuracy of the fitting.
c     do 210 i=1,nr
c     call cgntcs(roots(i),ew,ez,tc,ng)
c     ce=(0d0,0d0)
c     oo 220 l=1,ng
c 220 ce=ce+cm(l,j)*tc(l)
c 210 write(*,'(1x,a,6f12.5)') 'compare',roots(i),y(i),ce
      endif
   20 continue
c
c     --- Copy cm to the same l but different m states when
c         the calculation is one without the spin-orbit coupling.
      if(.not. ls) then
      mxlcmp=sqrt(dble(jmxcmp)-5d-1)+1
      do 100 j=1,mxlcmp
      i0=(j-1)**2+1
      do 100 i=i0+1,j**2
  100 call equarr(cm(1,i0),cm(1,i),ng)
      endif
      deallocate(roots,tc,y,c,d,extr)
      end
