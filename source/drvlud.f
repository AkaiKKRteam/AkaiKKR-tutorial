      subroutine drvlud(a,b,n)
c----------------------------------------------------------------------
c     Given a real*8 matrix a, this returns its inverse matrix.
c     After calling of this program a is replaced by the lu
c     decomposition. The half pivoting is employed.
c     coded by H.Akai, 1 Dec. 1996, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (nmax=500)
      real*8 a(n,n),b(n,n)
      integer indx(nmax)
      if(n .gt. nmax) call errtrp(1,'drvlud','n too large')
c     do 10 i=1,n
c     do 20 j=1,n
c  20 b(i,j)=(0d0,0d0)
c  10 b(i,i)=(1d0,0d0)
      call rludcm(a,n,indx,d)
      do 30 j=1,n
   30 call rlubks(a,n,indx,b(1,j))
      a(1,1)=a(1,1)*d
      end
