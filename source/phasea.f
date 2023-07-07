      subroutine phasea(v,tchc,tchs,fcs,mxl,mxlcmp,elvl,tm,ew,ez,ng
     &                 ,dr,xr,meshr,rstr,wk,sra,msg,corlvl,config
     &                 ,opnc,ebtm)
c-----------------------------------------------------------------------
c     Construct the Tchebycheff expansion of the phase function.
c     coded by H.Akai, 1983, julich
c     modified by H. Akai, 25 Aug. 1999, Osaka
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 tchc(ng,mxl**2),tchs(ng,mxl**2),v(meshr),dr(meshr)
     &      ,xr(meshr),fcs(4,mxl**2),elvl(ng),tm(ng,ng)
     &      ,rstr(meshr,mxl**2,ng),wk(meshr),corlvl(18),config(18)
     &      ,ec(2)
      real*8,allocatable::sn(:,:),cn(:,:),tc1(:),tc2(:),tc(:),td(:)
     &      ,tchr(:,:,:),rj(:),dj(:),rn(:),dn(:)
c     complex*16 t
      integer lc(18)
c     complex*16 t
      logical sra,msg,chk,opnc,bisect,valenc
      data c/274.0720442d0/,chk/.false./, itrmx/20/, bisect/.false./
     &    ,wdth/1d-4/
c            1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
     &  ,lc/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0/
      allocate(sn(mxl,ng),cn(mxl,ng),tc1(ng),tc2(ng),tc(ng),td(ng)
     &        ,tchr(ng,mxl,2),rj(mxl),dj(mxl),rn(mxl),dn(mxl))
      pi=4d0*atan(1d0)
      mxj=mxl**2
      lmx=mxlcmp-1
      rtin=xr(meshr-1)
      call clrarr(tchc,ng*mxj)
      call clrarr(tchs,ng*mxj)
      call clrarr(tchr,2*ng*mxl)
      do 10 k=1,ng
      ev=elvl(k)
      if(sra) ev=ev+(ev/c)**2
      call mdbsl(ev*rtin**2,lmx,rj,dj,rn,dn)
      do 10 j=1,mxlcmp
      jj=j
      i=(j-1)**2+1
      call radial(elvl(k),jj,rstr(1,i,k),rxg,dxg,v,v,dr,xr,meshr,sra)
      sn(j,k)=rtin**j*(rtin*rj(j)*dxg-dj(j)*rxg)
      cn(j,k)=rtin**(1-j)*(rtin*rn(j)*dxg-dn(j)*rxg)
c     if(j .eq. 3) then
c     ea=ev**(jj-1)
c     t=-cn(j,k)-(0d0,1d0)*ea*csqrt(cmplx(ev,1d-3))*sn(j,k)
c     t=-cn(j,k)/sn(j,k)/ea-(0d0,1d0)*csqrt(cmplx(ev,1d-3))
c     write(*,'(1p,10e15.7)')elvl(k),cn(j,k),sn(j,k),t
c     write(*,'(1p,10e15.7)')elvl(k),t
c     write(*,'(1p,10e15.7)')elvl(k),rxg,dxg
c     endif
      do 10 l=1,ng
      tchr(l,j,1)=tchr(l,j,1)+rxg*tm(l,k)
   10 tchr(l,j,2)=tchr(l,j,2)+dxg*tm(l,k)
      e1=ew-ez
      e2=ew+ez
      call gntcs(e1,ew,ez,tc1,ng)
      call gntcs(e2,ew,ez,tc2,ng)
      do 20 j=1,mxlcmp
      jj=j
      i=(j-1)**2+1
      r1=0d0
      r2=0d0
      do 30 l=1,ng
      r1=r1+tchr(l,j,1)*tc1(l)
   30 r2=r2+tchr(l,j,1)*tc2(l)
c     r1=1d0
c     r2=1d0
      e0=e2
c     e0=e1
      if(r1*r2 .lt. 0d0) then
c     write(*,*)'r1*r2<0 for j=',j
c     if(.false.) then
      e0=ew
      de=ez
      do 40 itr=1,4
      call gntcs(e0,ew,ez,tc,ng)
      rr=0d0
      do 50 l=1,ng
   50 rr=rr+tchr(l,j,1)*tc(l)
      sgn=1d0
      if(rr*r2 .gt. 0d0) sgn=-1d0
      de=5d-1*de
   40 e0=e0+sgn*de
      itrf=0
      do 60 itr=1,itrmx
      call gntcds(e0,ew,ez,tc,td,ng)
      rr=0d0
      rd=0d0
      do 70 l=1,ng
      rr=rr+tchr(l,j,1)*tc(l)
   70 rd=rd+tchr(l,j,1)*td(l)
      dlt=-rr/rd
      if(abs(dlt) .lt. 1d-5) exit
      e0=e0+dlt
   60 itrf=itr
      if(itrf .eq. itrmx) then
      call errtrp(2,'phasea','e0 not found')
      write(6,'(a,i2)')' l=',j-1
      e0=e2
      endif
      if(e0 .lt. e1 .or. e0 .gt. e2) then
      call errtrp(2,'phasea','illegal e0, converted')
      e0=e2
      endif
      if(msg .and. chk) write(*,'(a,i2,a,f12.7)')' j=',j,'  e0=',e0
c     write(*,'(a,i2,a,f12.7)')' j=',j,'  e0=',e0
      endif
c
      ev=e0
      if(sra) ev=ev+(ev/c)**2
      call mdbsl(ev*rtin**2,jj-1,rj,dj,rn,dn)
      call radial(e0,jj,wk,rxg,dxg,v,v,dr,xr,meshr,sra)
      fcs(1,i)=rtin**(1-j)*(rtin*rn(j)*dxg-dn(j)*rxg)
      fcs(2,i)=rtin**j*(rtin*rj(j)*dxg-dj(j)*rxg)
      fcs(3,i)=e0
      fcs(4,i)=0d0
      do 80 k=1,ng
      cn(j,k)=(cn(j,k)-fcs(1,i))/(elvl(k)-fcs(3,i))
      sn(j,k)=(sn(j,k)-fcs(2,i))/(elvl(k)-fcs(3,i))
      do 80 l=1,ng
      tchc(l,i)=tchc(l,i)+cn(j,k)*tm(l,k)
   80 tchs(l,i)=tchs(l,i)+sn(j,k)*tm(l,k)
      do 20 icor=18,1,-1
      if(config(icor) .gt. 0d0 .and. lc(icor)+1 .eq. jj
     &   .and. corlvl(icor) .lt. 0d0 .and. corlvl(icor) .gt. e1) then
c     --- such a core state should be treated as a partially
c          occupied core state.
      do 90 ic=1,2
      ec(ic)=corlvl(icor)
c     --- use Newton-Raphson to obtain the exact energy satisfying
c         R=0 and dR/dr=0, whose valued are given by Tchebychev
c         polinomials.
      itrf=0
      do 100 itr=1,itrmx
      call gntcds(ec(ic),ew,ez,tc,td,ng)
      rx=0d0
      dx=0d0
      do 110 n=1,ng
      rx=rx+tchr(n,jj,ic)*tc(n)
  110 dx=dx+tchr(n,jj,ic)*td(n)
      dec=-rx/dx
c     write(*,'(a,i3,4f12.7)')'itr,ec,rx,dec',itr,ec(i),rx,dec
      ec(ic)=ec(ic)+dec
      if(abs(dec) .lt. 1d-6) exit
  100 itrf=itr
      if(itrf .eq. itrmx) call errtrp(2,'phasea','ec not found')
c     --- ec(1) and ec(2) corresponding to the upper and lower bounds
c         of the band found.
c     write(*,'(a,f12.7)')'ec=',ec(i)
   90 continue
      ecor=corlvl(icor)
c
c     --- use a bisection to nawrrow the interval where the energy
c         satisfying c+is=0 exists. This pre-conditining may be used
c         if the Newton-Raphson meets any troubles.
c         exlwr is the value of ex at ec(2).
      if(bisect .and. ec(2) .lt. 0d0) then
      eupper=min(0d0,ec(1))
      elower=ec(2)
      call gntcs(elower,ew,ez,tc,ng)
      cx=0d0
      sx=0d0
      do 120 n=1,ng
      cx=cx+tchc(n,i)*tc(n)
  120 sx=sx+tchs(n,i)*tc(n)
      cx=cx*(ecor-e0)+fcs(1,i)
      sx=sx*(ecor-e0)+fcs(2,i)
      erel=ecor
      if(sra) erel=ecor+(ecor/c)**2
      eak=erel**(jj-1)*sqrt(-erel)
      exlwr=cx-sx*eak
      do 130 itr=1,10
      call gntcs(ecor,ew,ez,tc,ng)
      cx=0d0
      sx=0d0
      do 140 n=1,ng
      cx=cx+tchc(n,i)*tc(n)
  140 sx=sx+tchs(n,i)*tc(n)
      cx=cx*(ecor-e0)+fcs(1,i)
      sx=sx*(ecor-e0)+fcs(2,i)
      erel=ecor
      if(sra) erel=ecor+(ecor/c)**2
      eak=erel**(jj-1)*sqrt(-erel)
      ex=cx-sx*eak
      if(ex*exlwr .gt. 0d0) then
      elower=ecor
      else
      eupper=ecor
      endif
      ecor=5d-1*(elower+eupper)
      write(*,'(4f12.7)')ecor,ex
  130 continue
      endif
c
c     --- use Newton-Raphson to obtain the exact energy satisfying
c         c+is=0, where c and s are the valued given by Tchebychev
c         polinomials.
      itrf=0
      do 150 itr=1,itrmx
      valenc=ecor .gt. 0d0
      if(valenc) exit
      call gntcds(ecor,ew,ez,tc,td,ng)
      cx=0d0
      sx=0d0
      cxd=0d0
      sxd=0d0
      do 160 n=1,ng
      cx=cx+tchc(n,i)*tc(n)
      sx=sx+tchs(n,i)*tc(n)
      cxd=cxd+tchc(n,i)*td(n)
  160 sxd=sxd+tchs(n,i)*td(n)
c     --- be careful that cx and sx used for cxd and sxd
c         should be those before renormalization. For
c         this reason the calculation of cxd and sxd
c         comes first.
      cxd=cxd*(ecor-e0)+cx
      sxd=sxd*(ecor-e0)+sx
      cx=cx*(ecor-e0)+fcs(1,i)
      sx=sx*(ecor-e0)+fcs(2,i)
      erel=ecor
      if(sra) erel=ecor+(ecor/c)**2
      eak=erel**(jj-1)*sqrt(-erel)
      ex=cx-sx*eak
      exd=cxd-sxd*eak-sx*(dble(jj)-5d-1)*eak/erel
      dec=-ex/exd
c     write(*,'(a,i3,5f12.7)')'itr,ec,ex,exd,dec',itr,ecor,ex,exd,dec
      ecor=ecor+dec
      if(abs(dec) .lt. 1d-6) exit
  150 itrf=itr
      if(valenc) exit
      if(itrf .eq. itrmx) call errtrp(2,'phasea','ecor not found')
      corlvl(icor)=ecor
c     --- ecor satisfying c+is=0 found.
c     ecor=5d-1*(ec(1)+ec(2))
c     ecor=corlvl(icor)
c     wdth=5d-1*(ec(1)-ec(2))
c     wdth=1.2d0*wdth
c     wdth=0.3d0*wdth
c     --- wdth=0 corresponds to the procedure without smooth crossover.
c         Although there is no reason, wdth=1d-4 seems to work rather
c         well. This however might depend on the system.
      wt=-imag(log(dcmplx(ebtm-ecor,wdth)))/pi+1d0
      if((ecor. lt. 0d0) .or.
     &   (.not. opnc) .and. jj .gt. mxlcmp) then
      config(icor)=wt*dble(2*jj-1)
c     write(*,'(4f12.7)')ebtm,ecor,wdth,wt
c     write(*,'(a,5f12.7)')
c    & 'ecor,e2,e1,ebtm,wt',ecor,ec(2),ec(1),ebtm,wt
      exit
      endif
      endif
   20 continue
      do 170 j=1,mxlcmp
      i0=(j-1)**2+1
      do 170 i=i0+1,j**2
      do 180 n=1,ng
  180 call equarr(rstr(1,i0,n),rstr(1,i,n),meshr)
      call equarr(tchc(1,i0),tchc(1,i),ng)
      call equarr(tchs(1,i0),tchs(1,i),ng)
  170 call equarr(fcs(1,i0),fcs(1,i),4)
c     write(*,*)
c     do 180 n=1,ng
c     call gntcs(elvl(n),ew,ez,tc,ng)
c     ev=elvl(n)
c     if(sra) ev=ev+(ev/c)**2
c     eef=ev-fcs(3,5)
c     cx=0d0
c     sx=0d0
c     do 190 m=1,ng
c     cx=cx+tchc(m,5)*tc(m)
c 190 sx=sx+tchs(m,5)*tc(m)
c     cx=cx*eef+fcs(1,5)
c     sx=sx*eef+fcs(2,5)
c     t=-cx-(0d0,1d0)*csqrt(cmplx(ev,1d-3))*sx
c 180 write(*,'(1p,10e15.7)')elvl(n),t
c 180 write(*,'(1p,10e15.7)')elvl(n),cx,sx
c     write(*,*)
      deallocate(sn,cn,tc1,tc2,tc,td,tchr,rj,dj,rn,dn)
c     stop
      end
