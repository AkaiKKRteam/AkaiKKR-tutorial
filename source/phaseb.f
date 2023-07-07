      subroutine phaseb(v,tchc,tchs,fcs,mxl,mxlcmp,elvl,tm,ew,ez,ng,dr
     &                 ,xr,meshr,rstr,wk,sra,is,msg,corlvl,config,opnc
     &                 ,ebtm)
c-----------------------------------------------------------------------
c     -----------------------------------
c     --- spin-orbit included version ---
c     -----------------------------------
c     Construct the Tchebycheff expansion of the phase function.
c     coded by H.Akai, 1983, julich
c     This version includes the calculation of the spin-orbit coupling
c     constants.
c     coded by H.Akai, 1994, Osaka
c     adopted to spin-orbit version, Nov 1995, Osaka
c     modified by H. Akai, 25 Aug. 1999, Osaka
c     Modified by H. Akai, Tokyo, Jan. 24, 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v(meshr),dr(meshr),xr(meshr),wk(meshr,2*mxl+1)
     &      ,elvl(ng),tm(ng,ng),fcs(4,mxl**2),tchc(ng,mxl**2)
     &      ,tchs(ng,mxl**2),rstr(meshr,mxl**2,ng),corlvl(18),config(18)
      real*8,allocatable::sn(:,:),cn(:,:),tchr(:,:)
     &      ,tc1(:),tc2(:),tc(:),td(:),rj(:),dj(:),rn(:),dn(:)
      integer lc(18)
      logical msg,chk,test,sra,opnc,ncor
      data c/274.0720442d0/,chk/.false./,test/.false./, itrmx/100/
     &    ,wdth/1d-4/
c            1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
     &  ,lc/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0/
      pi=4d0*atan(1d0)
      mxj=mxl**2
      allocate(sn(mxj,ng),cn(mxj,ng),tchr(ng,mxj),tc1(ng),tc2(ng),tc(ng)
     &        ,td(ng),rj(mxl),dj(mxl),rn(mxl),dn(mxl))
c     if(i .eq. 2) then
c     test=.true.
c     else
c     test=.false.
c     endif
      if(test) call errtrp(3,'phaseb','spin-orbit suppressed')
      spn=(-1d0)**(is-1)
      lmax=mxlcmp-1
      mxjcmp=mxlcmp**2
      rtin=xr(meshr-1)
      v0=v(meshr)-v(meshr-1)
      if(chk) write(*,*)
      if(test) v0=0d0
      call clrarr(tchc,ng*mxj)
      call clrarr(tchs,ng*mxj)
      call clrarr(tchr,ng*mxj)
      do 10 k=1,meshr-1
   10 wk(k,1)=v(k)*xr(k)
      call diffn(wk(1,1),wk(1,2),dr,meshr-1)
      do 20 k=1,meshr-1
   20 wk(k,2)=(wk(k,2)-v(k))/xr(k)**2
      do 30 k=1,ng
      ev=elvl(k)
      if(sra) ev=ev+(ev/c)**2
      call mdbsl(ev*rtin**2,lmax,rj,dj,rn,dn)
      do 40 m=-lmax,lmax
      mm=3+lmax+m
      wk(meshr,mm)=v(meshr)
      fct=spn*dble(m)/c**2
c     --- spin moment 1/2 and factro 2 appearing in spin-orbit
c         coupling cancels to each other, remaining no extra
c         factor.
c     fct=spn*dble(-m)/c**2
c     write(*,*)k,m,fct
      if(test) fct=0d0
c     if(abs(m) .eq. 2) fct=0d0
      if(sra) then
      do 50 kk=1,meshr-1
   50 wk(kk,mm)=v(kk)+fct*wk(kk,2)/(1d0-(v(kk)-elvl(k))/c**2)**2
      else
c     --- for 'nrlls' case, the potential is m dependent but not depend
c         on energy.
              do 60 kk=1,meshr-1
   60 wk(kk,mm)=v(kk)+fct*wk(kk,2)
      endif
   40 continue
      do 30 j=1,mxjcmp
      jj=sqrt(dble(j)-5d-1)+1
      m=j-1-jj*(jj-1)
      mm=3+lmax+m
      call radial(elvl(k),jj,rstr(1,j,k),rxg,dxg,v,wk(1,mm),dr,xr
     &           ,meshr,sra)
c
c-------------------------------------------------------
c     if(.false.)then
c     if(sra) then
c     do 120 kk=1,meshr-1
c 120 wk(kk,1)=wk(kk,2)*(rstr(kk,j,k)/xr(kk))**2
c    &        /(1d0-(v(kk)-elvl(k))/c**2)**2
c     else
c     do 122 kk=1,meshr-1
c 122 wk(kk,1)=wk(kk,2)*(rstr(kk,j,k)/xr(kk))**2
c     endif
c     sum=0d0
c     do 124 kk=1,meshr-3,2
c     sum=sum+2d0*(wk(kk,1)*xr(kk)**2*dr(kk)
c    &       +4d0*wk(kk+1,1)*xr(kk+1)**2*dr(kk+1)
c    &       +wk(kk+2,1)*xr(kk+2)**2*dr(kk+2))/3d0/c**2
c 124 write(*,'(3x,i4,1p5e15.6)')kk+1,xr(kk+1),rstr(kk+1,j,k)
c    &    ,wk(kk+1,2)*xr(KK+1)**2,wk(kk+1,mm)*xr(kk+1),sum
c     slprm=2d0*fintgr(wk,dr,xr,meshr)/c**2
cc    slprm is increasing function with respect to energy because
cc    the nodes of higher energy radial wave function come
cc    more inside.
c     write(*,'(1x,a,f10.5,2i3,f10.5)')'  e,l,m,sl=',
c    &     elvl(k),jj,m,slprm
c     endif
c-------------------------------------------------------
c
      dxg=dxg+spn*dble(m)*v0*rxg/c**2
      sn(j,k)=rtin**jj*(rtin*rj(jj)*dxg-dj(jj)*rxg)
      cn(j,k)=rtin**(1-jj)*(rtin*rn(jj)*dxg-dn(jj)*rxg)
      do 30 l=1,ng
   30 tchr(l,j)=tchr(l,j)+rxg*tm(l,k)
      e1=ew-ez
      e2=ew+ez
      call gntcs(e1,ew,ez,tc1,ng)
      call gntcs(e2,ew,ez,tc2,ng)
      do 70 j=1,mxjcmp
      jj=sqrt(dble(j)-5d-1)+1
      m=j-1-jj*(jj-1)
      mm=3+lmax+m
      r1=0d0
      r2=0d0
      do 80 l=1,ng
      r1=r1+tchr(l,j)*tc1(l)
   80 r2=r2+tchr(l,j)*tc2(l)
      e0=e2
      if(r1*r2 .lt. 0d0) then
      e0=ew
      de=ez
      do 90 itr=1,4
      call gntcs(e0,ew,ez,tc,ng)
      rr=0d0
      do 100 l=1,ng
  100 rr=rr+tchr(l,j)*tc(l)
      sgn=1d0
      if(rr*r2 .gt. 0d0) sgn=-1d0
      de=5d-1*de
   90 e0=e0+sgn*de
      itrf=0
      do 110 itr=1,itrmx
      call gntcds(e0,ew,ez,tc,td,ng)
      rr=0d0
      rd=0d0
      do 120 l=1,ng
      rr=rr+tchr(l,j)*tc(l)
  120 rd=rd+tchr(l,j)*td(l)
      dlt=-rr/rd
      if(abs(dlt) .lt. 1d-5) exit
      e0=e0+dlt
  110 itrf=itr
      if(itrf .eq. itrmx) then
      call errtrp(2,'phaseb','e0 not found')
      write(*,'(a,i2)') '      for j=',j
      e0=e2
      endif
      if(e0 .lt. e1 .or. e0 .gt. e2) then
      call errtrp(2,'phaseb','illegal e0, converted')
      e0=e2
      endif
      if(msg .and. chk) then
      call errtrp(3,'phaseb',' ')
      write(*,'(6x,a,i1,a,f12.7)')'j=',j,'  e0=',e0
      endif
      endif
c
      ev=e0
      if(sra) then
      ev=ev+(ev/c)**2
      wk(meshr,mm)=v(meshr)
      fct=spn*dble(m)/c**2
      if(test) fct=0d0
c     if(abs(m) .eq. 2) fct=0d0
      do 130 kk=1,meshr-1
  130 wk(kk,mm)=v(kk)+fct*wk(kk,2)/(1d0-(v(kk)-e0)/c**2)**2
      endif
      call mdbsl(ev*rtin**2,jj-1,rj,dj,rn,dn)
      call radial(e0,jj,wk,rxg,dxg,v,wk(1,mm),dr,xr,meshr,sra)
      dxg=dxg+spn*dble(m)*v0*rxg/c**2
      fcs(1,j)=rtin**(1-jj)*(rtin*rn(jj)*dxg-dn(jj)*rxg)
      fcs(2,j)=rtin**jj*(rtin*rj(jj)*dxg-dj(jj)*rxg)
      fcs(3,j)=e0
      fcs(4,j)=0d0
c
c-------------------------------------------------------
      if(chk)then
      if(sra) then
      do 140 kk=1,meshr-1
  140 wk(kk,mm)=wk(kk,2)*(wk(kk,1)/xr(kk))**2
     &        /(1d0-(v(kk)-e0)/c**2)**2
      else
      do 150 kk=1,meshr-1
  150 wk(kk,mm)=wk(kk,2)*(wk(kk,1)/xr(kk))**2
      endif
      slprm=2d0*fintgr(wk(1,mm),dr,xr,meshr)/c**2
      write(*,'(1x,a,f10.5,2i3,f10.5)')'  e,l,m,sl=',
     &     e0,jj-1,m,slprm
      endif
c-------------------------------------------------------
      do 160 k=1,ng
      cn(j,k)=(cn(j,k)-fcs(1,j))/(elvl(k)-fcs(3,j))
      sn(j,k)=(sn(j,k)-fcs(2,j))/(elvl(k)-fcs(3,j))
      do 160 l=1,ng
      tchc(l,j)=tchc(l,j)+cn(j,k)*tm(l,k)
  160 tchs(l,j)=tchs(l,j)+sn(j,k)*tm(l,k)
c
      do 70 icor=18,1,-1
      if(config(icor) .gt. 0d0 .and. lc(icor)+1 .eq. jj
     &   .and. corlvl(icor) .lt. 0d0 .and. corlvl(icor) .gt. e1) then
      if(j .eq. (jj-1)**2+1) config(icor)=0d0
c     write(*,'(3i3,f12.7)')j,jj,icor,config(icor)
c     --- such a core state should be taken care of so as not to be
c         fully included as a core state.
      ecor=corlvl(icor)
      itrf=0
      do 170 itr=1,itrmx
      ncor=ecor .gt. 0d0
      if(ncor) exit
      call gntcds(ecor,ew,ez,tc,td,ng)
      cx=0d0
      sx=0d0
      cxd=0d0
      sxd=0d0
      do 180 n=1,ng
      cx=cx+tchc(n,j)*tc(n)
      sx=sx+tchs(n,j)*tc(n)
      cxd=cxd+tchc(n,j)*td(n)
  180 sxd=sxd+tchs(n,j)*td(n)
c     --- be careful that cx and sx used for cxd and sxd
c         should be those before renormalization. For
c         this reason the calculation of cxd and sxd
c         comes first.
      cxd=cxd*(ecor-e0)+cx
      sxd=sxd*(ecor-e0)+sx
      cx=cx*(ecor-e0)+fcs(1,j)
      sx=sx*(ecor-e0)+fcs(2,j)
      erel=ecor
      if(sra) erel=ecor+(ecor/c)**2
      eak=erel**(jj-1)*sqrt(-erel)
      ex=cx-sx*eak
      exd=cxd-sxd*eak-sx*(dble(jj)-5d-1)*eak/erel
      dec=-ex/exd
c     write(*,'(a,i3,5f12.7)')'itr,ec,ex,exd,dec',itr,ecor,ex,exd,dec
      ecor=ecor+dec
      if(abs(dec) .lt. 1d-6) exit
  170 itrf=itr
      if(ncor) exit
      if(itrf .eq. itrmx) call errtrp(2,'phaseb','ecor not found')
      corlvl(icor)=ecor
c     --- ecor satisfying c+is=0 found.
c     --- wdth=0 corresponds to the procedure without smooth crossover.
      wt=-imag(log(dcmplx(ebtm-ecor,wdth)))/pi+1d0
      if((ecor. lt. 0d0) .or.
     &   (.not. opnc) .and. jj .gt. mxlcmp) then
      config(icor)=config(icor)+wt
c     write(*,'(a,5f12.7)')
c    &    'ecor,ebtm,wt,config=',ecor,ebtm,wt,config(icor)
      exit
      endif
      endif
   70 continue
c     write(*,*)'phaseb'
c     write(*,'(1x,a,16f10.5)')' fcs(1)=',(fcs(1,j),j=1,mxlcmp**2)
c     write(*,'(1x,a,16f10.5)')' fcs(2)=',(fcs(2,j),j=1,mxlcmp**2)
c     write(*,'(1x,a,16f10.5)')' fcs(3)=',(fcs(3,j),j=1,mxlcmp**2)
c     write(*,'(1x,a,16f10.5)')' fcs(4)=',(fcs(4,j),j=1,mxlcmp**2)
c     write(*,'(1x,a,(5f10.5))')' cn=',(cn(5,k),k=1,ng)
c     write(*,'(1x,a,(5f10.5))')' sn=',(sn(5,k),k=1,ng)
c     write(*,'(1x,a,(5f10.5))')' tchc=',(tchc(k,5),k=1,ng)
c     write(*,'(1x,a,(5f10.5))')' tchc=',(tchs(k,5),k=1,ng)
      deallocate(sn,cn,tchr,tc1,tc2,tc,td,rj,dj,rn,dn)
      end
