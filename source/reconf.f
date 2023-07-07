      subroutine reconf(e,ebtm,config,rcnfg,tchc,tchs,fcs,ew,ez,ng
     &                 ,mxl,mxlcmp,jj,sra,ls,opnc,rearth)
c-----------------------------------------------------------------------
c     Given parameters that fit S and C functions, this program
c     retunrs a suitbale value for the occupation number 'config' of the
c     core state of an energy close to e.
c     Coded by H. Akai, Tokyo, Oct. 30, 2020`
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 tchc(ng,mxl**2),tchs(ng,mxl**2),fcs(4,mxl**2)
      real*8,allocatable::tc(:),td(:)
      logical sra,ls,opnc,valenc,bisect,rearth,supprs
c     --- a non-zero value of smear causes non-integer core
c         occupation numbers even for ASA cases. This does
c         not affect the results much but in the cases
c         that this may be a problem, set smear=0d0.
c         "supprs" controls whether the reconfiguration of the core
c         states should be taken or not.
      data c/274.0720442d0/, itrmx/20/,bisect/.true./
     &    ,smear/3d-3/, supprs/.true./
c    &    ,smear/1d-10/, supprs/.false./
c    &    ,smear/0d0/, suprrs/.false./
      pi=4d0*atan(1d0)
      e1=ew-ez
      e2=ew+ez
      if(opnc .and. rearth .and. jj .eq. 4) then
      rcnfg=abs(config)
      return
      else if(e .lt. e1) then
      rcnfg=abs(config)
      return
c     else if(e .gt. e2) then
      else if(e .gt. 0d0) then
      rcnfg=0d0
      return
      endif
c
      if(supprs) then
c     --- to suppress the procedure that reconfigure the occupations
c         of semi core states, activate the following 7 lines.
      if(e .lt. ebtm) then
      rcnfg=abs(config)
      return
      else
      rcnfg=0d0
      return
      endif
      endif
c
      allocate(tc(ng),td(ng))
      do 10 j=(jj-1)**2+1,jj**2
      if((ls .or. j .eq. (jj-1)**2+1) .and. (abs(config) .gt. 0d0)) then
      if(j .eq. (jj-1)**2+1) rcnfg=0d0
c     write(*,'(a,i3,3f12.5)')'j,e1,ecor,e2=',j,e1,ecor,e2
      ecor=e
c
c     --- use a bisection to narrow the interval where the energy
c         satisfying c+is=0 exists. This pre-conditioning may be used
c         if the Newton-Raphson meets any troubles.
c         exlwr is the value of ex at e1=ew-ez.
      if(bisect) then
      eupper=min(0d0,e2)
      elower=e1
      call gntcs(elower,ew,ez,tc,ng)
      cx=0d0
      sx=0d0
      do 20 n=1,ng
      cx=cx+tchc(n,j)*tc(n)
   20 sx=sx+tchs(n,j)*tc(n)
      eef=elower-fcs(3,j)
      if(abs(fcs(4,j)) .lt. 1d-6) then
      cx=cx*eef+fcs(1,j)
      sx=sx*eef+fcs(2,j)
      else
      gam=fcs(4,j)/(eef**2+fcs(4,j)**2)
      cx=cx*eef*gam+fcs(1,j)
      sx=sx*eef*gam+fcs(2,j)
      endif
      erel=elower
      if(sra) erel=elower+(elower/c)**2
      eak=erel**(jj-1)*sqrt(-erel)
      exlwr=cx-sx*eak
      do 30 itr=1,4
      call gntcs(ecor,ew,ez,tc,ng)
      cx=0d0
      sx=0d0
      do 40 n=1,ng
      cx=cx+tchc(n,j)*tc(n)
   40 sx=sx+tchs(n,j)*tc(n)
      eef=ecor-fcs(3,j)
      if(abs(fcs(4,j)) .lt. 1d-6) then
      cx=cx*eef+fcs(1,j)
      sx=sx*eef+fcs(2,j)
      else
      gam=fcs(4,j)/(eef**2+fcs(4,j)**2)
      cx=cx*eef*gam+fcs(1,j)
      sx=sx*eef*gam+fcs(2,j)
      endif
      erel=ecor
      if(sra) erel=ecor+(ecor/c)**2
      eak=erel**(jj-1)*sqrt(-erel)
      ex=cx-sx*eak
      if(ex*exlwr .gt. 0d0) then
      if(itr .eq. 1 .and. abs(ex) .gt. abs(exlwr)) then
c     write(*,'(a,2f12.5)')' ex, exlwr',ex,exlwr
      rcnfg=abs(config)
      return
      endif
      elower=ecor
      else
      eupper=ecor
      endif
      ecor=5d-1*(elower+eupper)
c     write(*,'(a,5f12.7)')
c    & 'elower,ecor,eupper,ex,exlwr=',elower,ecor,eupper,ex,exlwr
   30 continue
      endif
c     --- end bisection
c
      do 50 itr=1,itrmx
c     write(*,'(a,i3,f12.7)')'itr,ecor=',itr,ecor
      call gntcds(ecor,ew,ez,tc,td,ng)
      cx=0d0
      sx=0d0
      cxd=0d0
      sxd=0d0
      do 60 n=1,ng
      cx=cx+tchc(n,j)*tc(n)
      sx=sx+tchs(n,j)*tc(n)
      cxd=cxd+tchc(n,j)*td(n)
   60 sxd=sxd+tchs(n,j)*td(n)
c     --- be careful: cx and sx used for cxd and sxd
c         should be those before renormalization. For
c         this reason the calculation of cxd and sxd
c         comes first.
      eef=ecor-fcs(3,j)
      if(abs(fcs(4,j)) .lt. 1d-6) then
      cxd=cxd*eef+cx
      sxd=sxd*eef+sx
      cx=cx*eef+fcs(1,j)
      sx=sx*eef+fcs(2,j)
      else
      gam=fcs(4,j)/(eef**2+fcs(4,j)**2)
      cxd=(cxd*eef+cx*(2d0*gam*fcs(4,j)-1d0))*gam
      sxd=(sxd*eef+sx*(2d0*gam*fcs(4,j)-1d0))*gam
      cx=cx*eef*gam+fcs(1,j)
      sx=sx*eef*gam+fcs(2,j)
      endif
      erel=ecor
      if(sra) erel=ecor+(ecor/c)**2
      if(erel .gt. 0d0) then
c     write(*,'(a,2i3,f12.7)')'j,itr,erel=',j,itr,erel
      erel=-1d-3
c     stop
c     valenc=.true.
c     exit
      endif
      eak=erel**(jj-1)*sqrt(-erel)
      ex=cx-sx*eak
      exd=cxd-sxd*eak-sx*(dble(jj)-5d-1)*eak/erel
      dec=-ex/exd
c     write(*,'(a,2i3,4f12.7)')
c    &   'j,itr,ecor,ex,exd,dec',j,itr,ecor,ex,exd,dec
      ecor=ecor+dec
      if(ecor .gt. e2) ecor=e1
      if(ecor .lt. e1) ecor=-1d-2
      if(abs(dec) .lt. 1d-6) then
      valenc=ecor .gt. 0d0
      exit
      endif
   50 itrf=itr
c     If ecor > 0, it is not a core state, and hence, do nothiing.
      if(valenc) exit
c
      if(itrf .eq. itrmx) then
c     --- the procedure searching ecor fails.
      call errtrp(1,'reconf','ecor not found')
      endif
c
c     --- ecor satisfying c+is=0 found.
c     --- smear=0 corresponds to the procedure without smooth crossover.
      wt=-imag(log(dcmplx(ebtm-ecor,smear)))/pi+1d0
c     write(*,'(a,i3,5f12.7)')
c    &    'j,ecor,ebtm,wt=',j,ecor,ebtm,wt
      if(.not. ls) wt=wt*dble(2*jj-1)
      if((ecor. lt. 0d0) .or.
     &   (.not. opnc) .and. jj .gt. mxlcmp) then
      rcnfg=rcnfg+wt
      endif
      endif
   10 continue
      deallocate(tc,td)
      end
