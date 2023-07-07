      subroutine gsdatpnew(wk,a,vc,atmicp,r,anclr,corl,ro,rmt,dr,xr
     &                 ,iwk,itype,natm,ntyp,meshr,ncmp,ncmpx,conc)
c-----------------------------------------------------------------------
c    New version of 'gsdatp'. While the old version uses Matteiss'
c    prescription to construct the crystal version, this version
c    makes use of renormalize atom method, which is known to mimic
c    the atomic potential in the crystal fairly well.
c    Coded by H. Akai 19 Apr. 2009.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mxa=441,ncmpmx=100)
      real*8 wk(mxa,8),atmicp(3,natm),r(3,3),cnf(18)
     &      ,anclr(ncmpx),corl(18,ncmpx),conc(ncmpx)
     &      ,ro(meshr,ncmpx),rmt(ntyp)
     &      ,dr(meshr,ncmpx),xr(meshr,ncmpx)
      integer itype(natm),ncmp(ntyp)
      logical do(ncmpmx)
      data zero/1d-10/ , small/1d-7/
      pi=4d0*atan(1d0)
      if(ncmpx .gt. ncmpmx)
     &   call errtrp(1,'gsdatp','ncmpx too large')
c
c     ---define a radial mesh for each site and each component.
      do 10 i=1,ntyp
      do 10 j=1,ncmp(i)
      call jip(i,j,ji)
      do(ji)=.true.
      rtin=a*rmt(i)
c     ---'ratm' corresponds to the standard atomic volume.
      ratm=(3d0*qvolum(anclr(ji),2)/4d0/pi)**(1d0/3d0)
   10 call rmesha(1d-6,rtin,ratm,dr(1,ji),xr(1,ji),meshr)
      call clrarr(ro,meshr*ncmpx)
c     ---now construct charge densities.
      do 20 i=1,ntyp
      do 20 ic=1,ncmp(i)
      call jip(i,ic,ici)
      if(do(ici)) then
      call clrarr(cnf,18)
c     ---first perform atomic calculation.
      call atmicv(anclr(ici),cnf,corl(1,ici),wk(1,1),wk(1,2)
     &           ,wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7)
     &           ,wk(1,8),mxa,ier)
      do 30 j=i,ntyp
      do 30 jc=1,ncmp(j)
      call jip(j,jc,jcj)
      if(abs(anclr(jcj)-anclr(ici)) .lt. zero) then
      do(jcj)=.false.
      do 40 l=1,18
   40 corl(l,jcj)=corl(l,ici)
      do 50 k=1,meshr
   50 ro(k,jcj)=polint(xr(k,jcj),wk(1,3),wk(1,1),mxa)
      qins=fintgr(ro(1,jcj),dr(1,jcj),xr(1,jcj),meshr)
     &        +5d-1*dr(meshr,jcj)
     &            *(ro(meshr-1,jcj)*xr(meshr-1,jcj)**2
     &             +ro(meshr,jcj)*xr(meshr,jcj)**2)
      rnorm=anclr(jcj)/qins
      do 60 k=1,meshr
   60 ro(k,jcj)=ro(k,jcj)*rnorm
      endif
   30 continue
      endif
   20 continue
c     do 90 i=1,ncmpx
c     s1=fintgr(ro(1,i),dr(1,i),xr(1,i),meshr)
c  90 write(*,'(1x,a,i3,2f12.5)')
c    &  'component,Z,charge=',i,anclr(i),s1
      end
