      subroutine cinvrm(a,b,n,mx,m,ibsadr,kmx,convrg)
c----------------------------------------------------------------------
c     Given a complex matrix a, this returns its inverse matrix.
c     After calling of this program a is replaced by the lu
c     decomposition. No pivoting is employed.
c     In this version the operations running over the first index of
c     the matrix come up to the outermost loop, suitable for a big
c     matrix with rather small size regardint the first index.
c     coded by H.Akai., Feb. 4, 1992, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 a(n,n,mx),b(n,n,mx),dmy
      complex*16,allocatable::za(:,:),c(:),wk(:)
      integer,allocatable::ipiv(:)
      logical convrg(kmx)
      iwork=16*n
      allocate(za(n,n),c(n),ipiv(n),wk(iwork))
      pi=4d0*atan(1d0)
      do 10 k=1,m
      iadr=ibsadr+k
      kk=mod(iadr-1,kmx)+1
      if(.not. convrg(kk)) then
      call equarc(a(1,1,k),za,n**2)
      call zgetrf(n,n,za,n,ipiv,info)
      do 20 i=1,n
   20 c(i)=za(i,i)
      call zgetri(n,za,n,ipiv,wk,iwork,info)
      call equarc(za,b(1,1,k),n**2)
      ip=0d0
      do 40 i=2,n
   40 a(i,i,k)=c(i)/c(1)
      c(1)=(1d0,0d0)
      endif
   10 continue
      deallocate(za,c,wk,ipiv)
      end
