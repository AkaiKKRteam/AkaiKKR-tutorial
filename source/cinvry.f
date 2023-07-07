      subroutine cinvry(a,b,nx,n,m)
c----------------------------------------------------------------------
c     Given a complex matrix a, this returns its inverse matrix.
c     After calling of this program a is replaced by the lu
c     decomposition. No pivoting is employed.
c     In this version the operations running over the first index of
c     the matrix come up to the outermost loop, suitable for a big
c     matrix with rather small size regardint the first index.
c     coded by H.Akai, Feb. 4, 1992, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(nx,nx,m),b(nx,nx,m),c
      integer,allocatable::ipiv(:)
      allocate(ipiv(n))
      iwork=nx**2*m
      do 10 k=1,m
      call zgetrf(n,n,a(1,1,k),nx,ipiv,info)
   10 call zgetri(n,a(1,1,k),nx,ipiv,b,iwork,info)
      call equarc(a,b,nx*nx*m)
      deallocate(ipiv)
      end
