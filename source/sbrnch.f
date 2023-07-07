      subroutine sbrnch(a,m,k)
c-----------------------------------------------------------------------
c     Select a proper branch so as to keep continuity.
c     coded by H.Akai, 1986, Osaka
c     modified by H. Akai, Tokyo, 30 Oct. 2019
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(m,k),twopii
      pi=4d0*atan(1d0)
      twopi=2d0*pi
      twopii=cmplx(0d0,twopi)
      do 10 i=1,m
      a(i,1)=a(i,1)-twopii*dble(floor((dimag(a(i,1))+pi)/twopi))
      do 10 j=2,k
   10 a(i,j)=a(i,j)
     &      -twopii*dble(floor((dimag(a(i,j)-a(i,j-1))+pi)/twopi))
c     do 20 i=1,m
c     write(*,'(/1x,a)')'***sbrnch'
c  20 write(*,'(1x,i4,f12.5)')(j,-dimag(a(i,j))/pi,j=1,k)
      end
