      subroutine atmicv(anclr,cnf,e,ro,dr,xr,z0,z1,z2,wk,rc,msr,ier)
c-----------------------------------------------------------------------
c              ***** program header *****
c       This program calculates a self-consistent atomic
c       potential and constructes a crystal potential
c       from it according mattheiss' prescription.
c
c       coded by H.Akai on Sept.   1980
c       revised on Oct. 23,  1982  ( at Osaka)
c       revised on Aug.  3,  1983  ( at Juelich)
c       revised on Dec. 12,  1984  ( at Juelich)
c
c       new version on Jan. 10, 1986 (at Osaka)
c       Revised by H. Akai, Tokyo, Sep.17, 2020
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      real*8 cnf(18),e(18),dr(msr),xr(msr)
     &      ,z0(msr),z1(msr),z2(msr),wk(msr)
     &      ,rc(msr),ro(msr)
      integer npq(18),l(18),ncnf(18),iord(18)
      character*1 lsym(4)
      logical dsp,flg
      data dsp/.true./
c
c              1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
c     ------------------------------------------------------------------
c         j     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
c     ------------------------------------------------------------------
      data npq/ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 4, 5, 6, 6, 5, 6, 7/
     &     , l/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0/
     &   ,ncnf/ 2, 2, 6, 2, 6,10, 2, 6,10, 2, 6,14,10, 2, 6,14,10, 2/
c
c              1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d
     &   ,iord/ 1, 2, 3, 4, 5, 7, 6, 8,10, 9,11,14,12,13,15,18,16,17/
c
     &   ,lsym/'s','p','d','f'/
     &   , itrstp/30/ , tol/1d-6/, dmp/.4d0/
      call atmcor(anclr,anc,nc)
      if(anclr .lt. 1d-10) then
c     --- if the vacancy is specified, the atomic radial mesh and the
c         atomic energy level for a hydrogen atom are assumed for
c         computational convenience.
      call mshatm(ams,bms,dr,xr,1d0,msr,msr)
      call guesse(anclr,e)
      do 10 k=1,msr
      rc(k)=0d0
   10 ro(k)=0d0
      return
      endif
      anv=anclr
      flg=.true.
      do 20 i=1,nc
   20 if(cnf(i) .gt. 0d0) flg=.false.
      if(flg) then
      do 30 j=1,18
      i=iord(j)
      cnf(i)=min(anv,dble(ncnf(i)))
      anv=anv-dble(ncnf(i))
      if(anv .lt. 0d0) exit
   30 cnf(i)=ncnf(i)
      if(anv .gt. 0d0) call errtrp(1,'atmicv','nshell>18 not supported')
      endif
      nf=0
c     --- nf: last occupid shell
      do 40 j=1,18
   40 if(cnf(j) .gt. 0d0) nf=j
      call mshatm(ams,bms,dr,xr,anclr,msr,msr)
      call clrarr(z1,msr)
      call guessz(anclr,z1,xr,msr,msr)
      call guesse(anclr,e)
      call utimer(time,0)
      itrf=0
      do 50 itr=1,itrstp
      itr1=itr
      call cstatc(z1,rc,ro,e,nc,cnf,ams,bms,xr,msr,wk,ier,msr)
c     write(*,*)(e(i),i=1,nf)
      call poisnb(ro,z2,anclr,ams,bms,xr,msr,msr)
      do 60 k=msr,1,-1
      z2(k)=z2(k)-fldf(ro(k))*xr(k)
      if(z2(k) .lt. 2d0) kmerr=k
   60 z2(k)=max(2d0,z2(k))
      call eranlb(z0,z1,z2,rms,itr,dmp,dr,xr,msr,kmerr)
      if(abs(rms) .lt. tol) exit
   50 itrf=itr
      if(itrf .eq. itrstp) write(*,'(a)')'   not converge'
      if(dsp) then
      call utimer(time,itr1)
      write(*,'(a,f5.2)')'   nuclear charge=',anclr
      write(*,'(/t10,a,t18,a,t30,a/t5,35(''-''))')'nl','cnf','energy'
      do 70 i=1,nf
      if(abs(cnf(i)) .gt. 1d-10) then
      j=l(i)+1
      write(*,'(t10,i1,a1,t16,f6.3,4x,f12.4)')npq(i),lsym(j),cnf(i),e(i)
      endif
   70 continue
      write(*,'(//)')
      endif
      end
