      subroutine ludcmp(a,n,np,indx,d)
c-----------------------------------------------------------------------
c  (C) Copr. 1986-92 Numerical Recipes Software !+!).
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (tiny=1d-20)
      complex*16 a(np,np),sum,cdum
      real*8,allocatable::vv(:)
      integer indx(n)
      allocate(vv(n))
      d=1d0
      do 20 i=1,n
      aamax=0d0
      do 10 j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
   10 continue
c     if (aamax.lt.tiny) pause 'singular matrix in ludcmp'
      if (aamax.lt.tiny) call errtrp(1,'ludcmp','singular matrix')
      vv(i)=1d0/aamax
   20 continue
      do 90 j=1,n
      do 40 i=1,j-1
      sum=a(i,j)
      do 30 k=1,i-1
      sum=sum-a(i,k)*a(k,j)
   30 continue
      a(i,j)=sum
   40 continue
      aamax=0d0
      do 60 i=j,n
      sum=a(i,j)
      do 50 k=1,j-1
      sum=sum-a(i,k)*a(k,j)
   50 continue
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if (dum.ge.aamax) then
      imax=i
      aamax=dum
      endif
   60 continue
      if (j.ne.imax)then
      do 70 k=1,n
      cdum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=cdum
   70 continue
      d=-d
      vv(imax)=vv(j)
      endif
      indx(j)=imax
      if(abs(a(j,j)).lt.tiny) a(j,j)=tiny
      if(j.ne.n)then
      cdum=1d0/a(j,j)
      do 80 i=j+1,n
      a(i,j)=a(i,j)*cdum
   80 continue
      endif
   90 continue
      deallocate(vv)
      end
