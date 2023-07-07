      subroutine cemesh(ef,ewidth,edelt,ebtm,e,kmx)
c--------------------------------------------------------------------
c     Generate a semielliptic energy contour. Mesh points are located
c     following fermi's distribution function such that they are
c     distributed densely near the real axis.
c     The following are the examples of typical cases used in the past.
c       edelt/3d-4/, r/0.65d0 /, h/2d-1/, rc/4d-1/
c       edelt/1d-3/, r/0.70d0 /, h/2d-1/, rc/4d-1/
c       edelt/1d-3/, r/1.20d0 /, h/2d-1/, rc/4d-1/
c       edelt/3d-3/, r/1.40d0 /, h/2d-1/, rc/4d-1/
c       edelt/1d-4/, r/0.40d0 /, h/5d-1/, rc/4d-1/
c     coded by H.Akai, 1983, Juelich
c     latest version, 30 Nov.1997, Osaka
c     Scaling of the imaginary part is modified.
c     H. Akai, Tokyo, Dec. 6, 2022.
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 e(kmx)
      logical scaled
c     data h/5d-1/,rc/4d-1/,one/0.99999999d0/
c     data h/1.5d-1/,rc/4d-1/,one/0.99999999d0/
c     data h/1d-1/,rc/4d-1/,one/0.99999999d0/
c     data h/5d-3/,rc/4d-1/,one/0.99999999d0/
      data h/2d-1/,rc/4d-1/,one/0.99999999d0/
     &    ,scaled/.true./, rl/5d-1/
      if(kmx .le. 1) then
      e(1)=ef+(0d0,1d0)*edelt
      else
      r=ewidth/2d0
      if(scaled .and. r .gt. rl) then
      hs=h/r
      else
      hs=h/rl
      endif
      pi=4d0*atan(1d0)
      kc=rc*dble(kmx)+5d-1
      beta=log(pi/(edelt/hs/r)-1d0)/dble(kmx-kc)
      f=pi*(exp(beta*dble(1-kc))+1d0)*one
      do 10 k=1,kmx
      theta=f/(exp(beta*dble(k-kc))+1d0)
   10 e(k)=ef+r*dcmplx(cos(theta)-1d0,hs*sin(theta))
      endif
      ebtm=dble(e(1))
c     ebtm=ef-ewidth
c     de=ewidth/100d0
c     do 10 k=1,51
c  10 e(k)=ebtm+(edelt+dble(k-1)*1d-3)*(0d0,1d0)
c     do 20 k=52,151
c  20 e(k)=e(51)+de*dble(k-51)
c     do 30 k=152,201
c  30 e(k)=ef+(edelt+dble(201-k)*1d-3)*(0d0,1d0)
      return
      end
