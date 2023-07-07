      subroutine cinvrz(p,q,mx,m,n)
c----------------------------------------------------------------------
c     Given a complex matrix p, this returns its inverse matrix.
c     After calling of this program p is replaced by the lu
c     decomposition. Half pivoting is employed.
c     In this version the operations running over the first index of
c     the matrix come up to the outermost loop, suitable for a big
c     matrix with rather small size regardint the first index.
c     coded by H.Akai., Feb. 4, 1992, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 p(mx,mx,n),q(mx,mx,n)
      complex*16,allocatable::a(:,:),b(:,:)
      integer,allocatable::indx(:)
      allocate(a(m,m),b(m,m),indx(m))
      do 40 k=1,n
      do 10 i=1,m
      do 20 j=1,m
      a(i,j)=p(i,j,k)
   20 b(i,j)=(0d0,0d0)
   10 b(i,i)=(1d0,0d0)
      call ludcmp(a,m,m,indx,d)
      do 30 j=1,m
   30 call lubksb(a,m,m,indx,b(1,j))
      a(1,1)=a(1,1)*d
      do 40 i=1,m
      do 40 j=1,m
      p(i,j,k)=a(i,j)
   40 q(i,j,k)=b(i,j)
      deallocate(a,b,indx)
      end
