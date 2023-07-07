      subroutine diagrs(a,np,n,e,wrk)
c-----------------------------------------------------------------------
c     Diagonalization of a real symmetric n*n matrix 'a', for which
c     the upper half must be given.
c     In output, 'e' and 'a' contain eigenvalues and eigenvectors,
c     respectively; k-th column of 'a' returns the normalized
c     eigenvector corresponding to e(k). 'wrk' is a working area.
c     eigenvalues are sorted into ascending order.
c     Coded by H.Akai 29 Dec. 1990
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(np,n),e(n),wrk(n)
      call trdiag(a,np,n,e,wrk)
      call tqlalg(a,np,n,e,wrk)
      call eigsrt(a,np,n,e)
c     call dspdia(a,np,n,e)
      end
      subroutine trdiag(a,np,n,d,e)
c-----------------------------------------------------------------------
c     Householder reduction of 'a' real symmetric matrix 'a' into
c     tridiagonal matrix. 'd' and 'e' contain the diagonal and the
c     subdiagonal components of the tridiagonal matrix.
c     Upper half of the matrix must be given.
c     See numerical recipes by h. press et. al., p.349.
c     Coded by H.Akai 29 Dec. 1990
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(np,n),d(n),e(n)
      data small/1d-30/
      do 70 i=n,2,-1
      l=i-1
      h=0d0
      scale=0d0
      if(l .gt. 1) then
      do 10 k=1,l
   10 scale=scale+abs(a(k,i))
      if(abs(scale) .lt. small) then
      e(i)=a(l,i)
      else
      do 20 k=1,l
      a(k,i)=a(k,i)/scale
   20 h=h+a(k,i)**2
      f=a(l,i)
      g=-sign(sqrt(h),f)
      e(i)=scale*g
      h=h-f*g
      a(l,i)=f-g
      f=0d0
      do 50 j=1,l
      a(i,j)=a(j,i)/h
      g=0d0
      do 30 k=1,j
   30 g=g+a(k,j)*a(k,i)
      do 40 k=j+1,l
   40 g=g+a(j,k)*a(k,i)
      e(j)=g/h
   50 f=f+e(j)*a(j,i)
      hh=f/(h+h)
      do 60 j=1,l
      f=a(j,i)
      g=e(j)-hh*f
      e(j)=g
      do 60 k=1,j
   60 a(k,j)=a(k,j)-f*e(k)-g*a(k,i)
      endif
      else
      e(i)=a(l,i)
      endif
   70 d(i)=h
      d(1)=0d0
      e(1)=0d0
      do 100 i=1,n
      l=i-1
      if(abs(d(i)) .gt. 0d0) then
      do 90 j=1,l
      g=0d0
      do 80 k=1,l
   80 g=g+a(k,i)*a(k,j)
      do 90 k=1,l
   90 a(k,j)=a(k,j)-g*a(i,k)
      endif
      d(i)=a(i,i)
      a(i,i)=1d0
      do 100 j=1,l
      a(i,j)=0d0
  100 a(j,i)=0d0
      end
      subroutine tqlalg(a,np,n,d,e)
c-----------------------------------------------------------------------
c     QL algorithm with implicit shift
c     d: diagonal part
c     e: subdiagonal part
c     In output, 'd' and 'a' contain eigenvalues and eigenvectors.
c     See numerical recipes by h. press et. al., p.360.
c     Coded by H.Akai 29 Dec. 1990
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(np,n),d(n),e(n)
      data small/1d-32/
      do 10 i=2,n
   10 e(i-1)=e(i)
      e(n)=0d0
      do 60 l=1,n
      itr=0
   20 do 30 m=l,n-1
      dd=abs(d(m))+abs(d(m+1))
      de=abs(e(m))+dd
      if(abs(dd-de) .lt. small) go to 40
   30 continue
      m=n
   40 if(m .ne. l) then
      if(itr .eq. 30) then
      write(6,1000)
 1000 format('   ***err in tqlalg...too many iteration')
      stop
      endif
      itr=itr+1
      g=5d-1*(d(l+1)-d(l))/e(l)
      r=sqrt(g**2+1d0)
      g=d(m)-d(l)+e(l)/(g+sign(r,g))
      s=1d0
      c=1d0
      p=0d0
      do 50 i=m-1,l,-1
      f=s*e(i)
      b=c*e(i)
      if(abs(f) .ge. abs(g)) then
      c=g/f
      r=sqrt(c**2+1d0)
      e(i+1)=f*r
      s=1d0/r
      c=c*s
      else
      s=f/g
      r=sqrt(s**2+1d0)
      e(i+1)=g*r
      c=1d0/r
      s=s*c
      endif
      g=d(i+1)-p
      r=(d(i)-g)*s+2d0*c*b
      p=s*r
      d(i+1)=g+p
      g=c*r-b
      do 50 k=1,n
      f=a(k,i+1)
      a(k,i+1)=s*a(k,i)+c*f
   50 a(k,i)=c*a(k,i)-s*f
      d(l)=d(l)-p
      e(l)=g
      e(m)=0d0
      go to 20
      endif
   60 continue
      end
      subroutine eigsrt(a,np,n,d)
c-----------------------------------------------------------------------
c     Sorts the eigenvalues into ascending order.
c     See numerical recipes by h. press et. al., p.348.
c     Coded by H.Akai 29 Dec. 1990
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(np,n),d(n)
      do 30 i=1,n-1
      k=i
      p=d(i)
      do 10 j=i+1,n
      if(d(j) .lt. p) then
      k=j
      p=d(j)
      endif
   10 continue
      if(k .ne. i) then
      d(k)=d(i)
      d(i)=p
      do 20 j=1,n
      p=a(j,k)
      a(j,k)=a(j,i)
   20 a(j,i)=p
      endif
   30 continue
      end
      subroutine dspdia(a,np,n,d)
c-----------------------------------------------------------------------
c     Displays eigenvalues and corresponding eigenvectors.
c     Coded by H.Akai 29 Dec. 1990
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 a(np,n),d(n)
      do 10 i=1,n-3,4
      write(6,'(//5x,4(a10,5x))')('eigenvalue',k=i,i+3)
      write(6,'(2x,4(1pd15.7))')(d(k),k=i,i+3)
      write(6,'(/5x,4(a11,4x))')('eigenvector',k=i,i+3)
      do 10 k=1,n
   10 write(6,'(2x,4(1pd15.7))')(a(k,j),j=i,i+3)
      if(mod(n,4) .ne. 0) then
      write(6,'(//5x,4(a10,5x))')('eigenvalue',k=n/4*4+1,n)
      write(6,'(2x,4(1pd15.7))')(d(k),k=n/4*4+1,n)
      write(6,'(/5x,4(a11,4x))')('eigenvector',k=n/4*4+1,n)
      do 20 k=1,n
   20 write(6,'(2x,4(1pd15.7))')(a(k,j),j=n/4*4+1,n)
      endif
      end
