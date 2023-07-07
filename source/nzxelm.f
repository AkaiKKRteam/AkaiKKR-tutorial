      subroutine nzxelm(urotat,mxl,isymop,natm,irotat,lmxtyp,itype,iatmp
     &                 ,type,mtrx,iatm,ntyp,ncmpx,mtabl,nznum)
c-----------------------------------------------------------------------
c     This program generate the table of non-zero element of
c     G(l,l',itype) that are used in the conductivity calculations.
c     Here, 'itype' indicates the type of sites.
c     A version of 'nzelem' coded by H. Akai, 7 Dec. 2008
c     Coded by H. Akai, 11 June 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 urotat((2*mxl-1)**2,mxl,24),u
      complex*16,allocatable::a(:),b(:),c(:,:,:),table(:)
      real*8,allocatable::r(:,:)
      integer isymop(24,2),irotat(natm,24,2),lmxtyp(*),itype(natm)
     &       ,iatmp(natm),mtrx(mxl**2,mxl**2,ncmpx),iatm(ntyp)
     &       ,mtabl(3,mxl**4*ncmpx)
      character*4 buff(256),format*256,type(ntyp)*8
      data zero/1d-10/,iprint/1/
      mmxl=mxl**2
      msiz=0
      do 10 n=1,natm
      mmx=(lmxtyp(itype(n))+1)**2
   10 msiz=msiz+(lmxtyp(itype(n))+1)**2
      allocate(a(msiz),b(msiz),c(mmxl,mmxl,ntyp),table(mmxl*ntyp)
     &         ,r(msiz,2),stat=ierr)
      if(ierr .ne. 0) call errtrp(1,'nztelm','allocation error')
      call clrarc(a,msiz)
      call clrarc(b,msiz)
      call clrarc(c,mmxl**2*ntyp)
      call getnum(1011,r,2*msiz)
c     --- msiz is the size needed to cover all the atoms in the
c         unit cell.
c     --- fill in the array "a" with random numbers.
      do 20 i=1,msiz
   20 a(i)=dcmplx(r(i,1),r(i,2))
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
c     --- we also have to take account that the atom positions
c         are also rotated.
c     --- site i is occupied by atom irotat(i,iop) after the rotation.
      j=irotat(i,iop,ip)
c     iatmph=(iatmp(i)-1)/2
c     jatmph=(iatmp(j)-1)/2
      iatmph=iatmp(i)-1
      jatmph=iatmp(j)-1
c     write(*,*)iop,ip,i,j,iatmph,jatmph,mr,mc,p
   50 if(l .le. lmxtyp(itype(i))+1)
     &    b(iatmph+mc)=b(iatmph+mc)+a(jatmph+mr)*u
c     --- Now, b(i)=Sum_k a(k) * conjg(U(i,k))
      endif
   40 continue
      do 60 i=1,ntyp
      mxmxl=(lmxtyp(i)+1)**2
      n=iatm(i)
      ibs=iatmp(n)-1
      do 60 mc=1,mxmxl
      do 60 mr=1,mxmxl
   60 c(mr,mc,i)=c(mr,mc,i)+conjg(b(ibs+mr))*b(ibs+mc)
      endif
   30 continue
      nznum=0
      nzero=0
      do 80 ji=1,ncmpx
      call jipinv(ji,i,j)
      jlst=nznum
      mmx=(lmxtyp(i)+1)**2
      do 80 n1=1,mmx
      do 80 n2=n1,mmx
      mtrx(n1,n2,ji)=0
      mtrx(n2,n1,ji)=0
      if(abs(c(n1,n2,i)) .gt. zero) then
      nzero=nzero+1
      l=nznum+1
      do 90 k=jlst+1,nznum
      if(abs(c(n1,n2,i)-table(k)) .lt. zero) l=k
   90 if(abs(c(n1,n2,i)+table(k)) .lt. zero) l=-k
      if(l .eq. nznum+1) then
      nznum=nznum+1
      table(nznum)=c(n1,n2,i)
      mtabl(1,nznum)=n1
      mtabl(2,nznum)=n2
      mtabl(3,nznum)=ji
      endif
      mtrx(n1,n2,ji)=l
      mtrx(n2,n1,ji)=l
      endif
   80 continue
      
c     write(*,*)
c     write(*,'(1x,a,i5)')'table size=',nznum
c     write(*,'(1x,a,i5)')'number of nonzero elements=',nzero
c     write(22)mtrx
c     write(*,'((1x,1p,6e13.5))')(table(j)/dble(iw),j=1,nznum)
c
c     ------------- optional print-out ---------------
      if(iprint .eq. 1) then
c     write(*,'(1x,36i3)')((mtrx(i,j),i=1,msiz),j=1,msiz)
      do 100 ji=1,ncmpx
      call jipinv(ji,i,j)
      mmx=(lmxtyp(i)+1)**2
      write(format,'(''(5x,'',i3,''(''''-''''))'')')2+4*mmx
      write(*,'(/8x,3a,i3)')'type=',type(i),'  component=',j
      write(*,format)
      write(format,'(''(1x,'',i3,''a4)'')')mmx+1
      do 100 l=1,mmx
      do 110 k=1,l
  110 buff(k)=' '
      do 120 k=l,mmx
      buff(k+1)='   .'
  120 if(mtrx(l,k,ji) .ne. 0) write(buff(k+1),'(i4)')mtrx(l,k,ji)
  100 write(*,format)(buff(k),k=1,mmx+1)
      write(*,*)
c     write(*,'(4i3)')(i,(mtabl(j,i),j=1,3),i=1,nznum)
      endif
      deallocate(a,b,c,table,r,stat=ierr)
      if(ierr .ne. 0) call errtrp(1,'ntdelm','deallocation error')
c     stop
      end
