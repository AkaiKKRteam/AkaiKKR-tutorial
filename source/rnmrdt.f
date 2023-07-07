      function rnmrdt(a)
c----------------------------------------------------------------------
c     Given a character string, this program returns the first
c     numeric data contained in the string. If the data does not
c     contain any feasible numerical data, the program returns 1d30.
c     Coded by H. Akai, 8 June 2016, Tokyo
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*1 a*(*),x,s,t,u
      character*1,allocatable::b(:),c(:)
      logical e,d,p,n,np,nd,ne,c1,c2,c3
      e(x)=x .eq. 'e' .or. x .eq. 'd' .or. x .eq. 'E' .or. x .eq. 'D'
      n(x)=x .ge. '0' .and. x .le. '9'
      p(x)=x .eq. '+' .or. x .eq. '-'
      d(x)=x .eq. '.'
      np(x)=n(x) .or. p(x)
      nd(x)=n(x) .or. d(x)
      ne(x)=n(x) .or. e(x)
      m=len_trim(a)
      allocate(b(m+1),c(m+1))
      b=' '
      do 10 i=1,m
   10 b(i+1)=a(i:i)
      c=b
      do 20 i=2,m+1 
      s=b(i-1)
      t=b(i)
      u=b(i+1)
      c1=nd(s) .and. e(t) .and. np(u)
      c2=(np(s) .or. ne(u)) .and. d(t)
      c3=(e(s) .or. nd(u)) .and. p(t)
   20 if(.not. (n(t) .or. c1 .or. c2 .or. c3)) b(i)=' '
      rnmrdt=rnmrdx(b,m+1,ierr)
      if(ierr .eq. 0) then
      do 30 i=1,m+1
   30 if(b(i) .eq. c(i)) c(i)=' '
      do 40 i=1,m
   40 a(i:i)=c(i+1)
      else
      rnmrdt=1d30
      endif
      deallocate(b,c)
      end 
