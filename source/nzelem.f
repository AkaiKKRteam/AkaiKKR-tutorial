      subroutine nzelem(urotat,mxl,isymop,natm,irotat,lmxtyp,itype
     &                 ,iatmp,atmtyp)
c-----------------------------------------------------------------------
c     This program generate the table of non-zero element of G and Xi
c     that are used in the conductivity calculations.
c     Coded by H. Akai, 7 Dec. 2008
c     Make a table for the non-zero elements of the single site matrix.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 urotat((2*mxl-1)**2,mxl,24),u
      complex*16,allocatable::a(:),b(:),c(:,:),table(:)
      real*8,allocatable::r(:)
      integer isymop(24,2),irotat(natm,24,2),lmxtyp(*),itype(natm)
     &       ,iatmp(natm)
      integer,allocatable::mtrx(:,:)
      character*4 buff(256),format*256,atmtyp(natm)*8,title*256
      logical skip
      data zero/1d-8/,iprint/0/, skip/.true./
      if(skip) return
      msiz=0
      do 10 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
   10 msiz=msiz+mmx
      allocate(a(msiz),b(msiz),c(msiz,msiz),table(msiz**2),r(msiz)
     &        ,mtrx(msiz,msiz))
      call clrarc(a,msiz)
      call clrarc(b,msiz)
      call clrarc(c,msiz**2)
      call getnum(1011,r,msiz)
c     --- msiz is the size spanned by all the atoms in the unit cell.
c     --- fill in the array "a" with random numbers.
      do 20 i=1,msiz
   20 a(i)=dcmplx(r(i),0d0)
c     --- the matrix simulating a general matrix, its row and column
c         being (l,m) and atomic sites in the unit cell, is given by the
c         direct product a x a.
      iw=0
      do 30 ip=1,2
      do 30 iop=1,24
      if(isymop(iop,ip) .eq. 1) then
      iw=iw+1
c     --- only rotations compatible with the crystal symmetry
c         are considered.
c     --- we first calculate wk1=f*U^+, where U^+ is the Hermite
c         conjugate of the rotation matrix U.
      call clrarc(b,msiz)
      do 40 l=1,mxl
c     --- for the case of improper rotations
      p=(-1d0)**((ip-1)*(l-1))
      mx=2*l-1
      lb=(l-1)**2
      do 40 m1=1,mx
      m0=mx*(m1-1)
      do 40 m2=1,mx
      m2m1=m0+m2
      if(abs(urotat(m2m1,l,iop)) .gt. zero) then
c     --- only non-zero elements of the rotation matrix are
c         taken into account. In order to handle this procedure
c         efficiently, the order of summation appraring in
c         the matrix product is changed from a straightforward
c         way.
c     --- take account of the effects of the inversion.
      u=p*urotat(m2m1,l,iop)
      u=conjg(u)
      mr=lb+m1
      mc=lb+m2
      do 50 i=1,natm
c     --- we also have to take account of that the atom positions
c         are also rotated.
c     --- site i is occupied by atom irotat(i,iop) after the rotation.
      j=irotat(i,iop,ip)
c     iatmph=(iatmp(i)-1)/2
c     jatmph=(iatmp(j)-1)/2
      iatmph=iatmp(i)-1
      jatmph=iatmp(j)-1
c     write(*,*)iop,ip,i,j,iatmph,jatmph,mr,mc,p
   50 if(l .le. lmxtyp(itype(i))+1) 
     &    b(iatmph+mc)=
     &    b(iatmph+mc)+a(jatmph+mr)*u
c     --- Now, b(i)=Sum_k a(k) * conjg(U(i,k))
      endif
   40 continue
      do 60 i=1,msiz
      do 60 j=1,msiz
   60 c(i,j)=c(i,j)+conjg(b(i))*b(j)
      endif
   30 continue
c     write(*,'(1x,9f10.3)')((abs(c(i,j)),i=1,msiz),j=1,msiz)
      jj=1
      table(1)=c(1,1)
      nzero=0
      do 80 n1=1,msiz
      do 80 n2=n1,msiz
      mtrx(n1,n2)=0
      mtrx(n2,n1)=0
      if(abs(c(n1,n2)) .gt. 1d-6) then
      nzero=nzero+1
      i=jj+1
      do 90 k=1,jj
      if(abs(c(n1,n2)-table(k)) .lt. 1d-6) i=k
   90 if(abs(c(n1,n2)+table(k)) .lt. 1d-6) i=-k
      if(i .eq. jj+1) then
      jj=jj+1
      if(jj .gt. msiz**2) call errtrp(1,'nzelem','table full')
      table(jj)=c(n1,n2)
      endif
      mtrx(n1,n2)=i
      mtrx(n2,n1)=i
      endif
   80 continue
c     write(*,*)
c     write(*,'(1x,a,i5)')'table size=',jj
c     write(*,'(1x,a,i5)')'number of nonzero elements=',nzero
c     write(22)mtrx
c     write(*,'((1x,1p,6e13.5))')(table(j)/dble(iw),j=1,jj)
c
c     ------------- optional print-out ---------------
      if(iprint .eq. 1) then
c     write(*,'(1x,36i3)')((mtrx(i,j),i=1,msiz),j=1,msiz)
      title=' '
      ip=9
      do 70 i=1,natm
      title(ip:ip+8)=atmtyp(i)
   70 ip=ip+4*(lmxtyp(itype(i))+1)**2
      write(*,'(/a)')title(1:ip-4)
      do 130 i=7,ip-2
  130 title(i:i)='-'
      ip=9
      do 140 i=1,natm
      title(ip:ip)='+'
  140 ip=ip+4*(lmxtyp(itype(i))+1)**2
      write(*,'(a)')title(1:ip-2)
c     write(format,'(''(5x,'',i3,''(''''-''''))'')')2+4*msiz
c     write(*,format)
      write(format,'(''(1x,'',i3,''a4)'')')msiz+1
      do 100 i=1,msiz
      do 110 j=1,i
  110 buff(j)=' '
      do 120 j=i,msiz
      buff(j+1)='   .'
  120 if(mtrx(i,j) .ne. 0) write(buff(j+1),'(i4)')mtrx(i,j)
  100 write(*,format)(buff(j),j=1,msiz+1)
      write(*,*)
      endif
      deallocate(a,b,c,table,r,mtrx)
c     stop
      end
