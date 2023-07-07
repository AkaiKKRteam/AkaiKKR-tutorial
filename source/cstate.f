      subroutine cstate(v,ro,rorg,corlvl,anc,dr,xr,match
     &                 ,config,esic,meshr,wk,sra,sdftyp,ebtm,asa,mxlcmp
     &                 ,openc,mxl,tchc,tchs,fcs,ew,ez,ng,ls,rcnfg)
c----------------------------------------------------------------------
c     Calculate core state by matching boundary condition.
c     coded by H.Akai, 1983, Juelich
c     very minor modification for message output.
c     by H. Akai, Feb 2004.
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v(meshr),dr(meshr),xr(meshr),ro(meshr),wk(meshr,6)
     &         ,rorg(20),npq(18),corlvl(18),config(18),rcnfg(18)
     &         ,tchc(ng,mxl**2),tchs(ng,mxl**2),fcs(4,mxl**2)
      integer match(18),l(18)
      logical sic,sicon,init,asa,ifkey,start,start2,on(6),openc,rearth
     &       ,sra,opnc(6),ls
      character sdftyp*12,sicor*24,vsymbl(6)*1
      data amagic/1d0/, start/.true./, start2/.true./
     &    ,istop/50/, tol/1d-8/, eb/-20d0/, sic/.false./
     &    ,on/3*.false.,3*.true./
     &    ,vsymbl/'s','p','d','f','g','h'/
c
c            1s 2s 2p 3s 3p 3d 4s 4p 4d 5s 5p 4f 5d 6s 6p 5f 6d 7s
     &  ,npq/ 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 4, 5, 6, 6, 5, 6, 7/
     &    ,l/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0/
      save start,start2
c
      anclr=0d0
      do 10 k=1,meshr/8,5
   10 anclr=max(anclr,-v(k)*xr(k))
      nclr=5d-1*anclr+5d-1
      rearth=(nclr-57)*(nclr-71)*(nclr-89)*(nclr-103) .le. 0
      if(start .and. sic) then
      sicor=' '
      do 20 i=1,mxl
   20 if(on(i)) sicor=trim(sicor)//', '//vsymbl(i)
      if(sicor .ne. ' ') then
      sicor=sicor(1:index(sicor,',')-1)
     &    //sicor(index(sicor,',')+1:len(sicor))
      if(index(sicor,',') .gt. 0)
     &   sicor=sicor(1:index(sicor,',',back=.true.)-1)
     &   //' and'//sicor(index(sicor,',',back=.true.)+1:len(sicor))
      write(*,'(/3x,a)')'***************************'
      write(*,'(5x,2a)')'core SIC for',sicor(1:len(sicor))
      write(*,'(5x,a,f6.3)')'amagic=',amagic
      write(*,'(3x,a/)')'***************************'
      endif
      start=.false.
      endif
      if(start2 .and. openc .and. rearth) then
      write(*,'(/3x,a)')'***************************'
      write(*,'(8x,a)')'open f-core'
      write(*,'(3x,a/)')'***************************'
      start2=.false.
      endif
c
c     --- if sic is on, openc also must be on.
      do 30 i=1,6
   30 opnc(i)=openc .or. (sic .and. on(i))
      call clrarr(ro,meshr)
c     write(*,'((1x,1p6e13.6))')(v(k),k=1,meshr,10)
c     write(*,'((1x,1p6e13.6))')(-5d-1*v(k)*xr(k),k=1,60)
      anc=0d0
      tint=0d0
      esic=0d0
      do 40 j=1,18
      rcnfg(j)=0d0
      if(abs(config(j)) .gt. 0d0) then
      jj=l(j)+1
      node=npq(j)-jj
      e=corlvl(j)
c     ei=e
      eold=e
      init=match(j) .lt. 1
      do 50 k=1,meshr
   50 wk(k,3)=v(k)
      do 60 ip=1,istop
      emax=1d10
      emin=-1d10
      if(e .gt. 1d2) e=0.1d0
c     if(abs(e) .lt. 1d-30) e=-0.1d0
      do 70 itr=1,istop
c     if(j .eq. 12)
c    &write(*,'(1x,a,i2,f12.5,2i4)')'itr,e,match,nn=',itr,e,match(j),nn
      call corada(e,jj,wk(1,1),rin,match(j),g1,g2,nn,wk(1,3)
     &           ,dr,xr,meshr,sra,asa)
c     write(*,'(1x,a,i2,f12.5,2i4,1p,2e13.6,i4)')
c    &  'itr,e,match,nn=',itr,e,match(j),nn
c    &         ,emin-eb,emax-eb,node
      if(nn .gt. node) then
      emax=min(e+eb,emax)
      e=max(emax*1.25d0,(emax+emin)*5d-1)-eb
      if(init) match(j)=0
      go to 70
      endif
      if(nn .lt. node) then
      emin=max(e+eb,emin)
      e=min(emin*0.75d0,(emax+emin)*5d-1)-eb
      if(init) match(j)=0
      go to 70
      endif
      dlt=-wk(match(j),1)*(g1-g2)
c     if(j .eq. 1)
c    & write(*,'(a,i4,1p5e20.12)')'match,wk,g1,g2,dlt,e',
c    & match(j),wk(match(j),1),g1,g2,dlt,e
      if(abs(dlt) .lt. tol) go to 80
      if(dlt .gt. 0d0) emin=max(emin,e+eb)
      if(dlt .lt. 0d0) emax=min(emax,e+eb)
      e=e+dlt
   70 continue
c     if(e .lt. 0d0) then
      if(e .lt. ebtm) then
      write(*,'(a,i3,a,i3)')
     & ' ***wrn in cstate...no convergence for nclr=',nclr,' j=',j
      write(*,'(1x,6f12.5)')(wk(k,1),k=1,meshr,10)
      else
      e=1d3
      endif
c
   80 do 90 k=1,meshr
      wk(k,6)=v(k)-wk(k,3)
      wk(k,1)=wk(k,1)/xr(k)**2
   90 wk(k,2)=wk(k,1)-1d-10
      sicon=sic .and.  on(jj)
      if(.not. sicon) go to 100
c     write(*,*)anclr,j,e
      call poisna(wk(1,1),wk(1,3),0d0,dr,xr,meshr)
      do 110 k=1,meshr
  110 wk(k,5)=5d-1*wk(k,3)
      if(ifkey('vbh',sdftyp)) then
      call excvbh(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('mjw',sdftyp)) then
      call excmjw(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('vwn',sdftyp)) then
      call excvwn(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('lmm',sdftyp)) then
      call exclmm(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('pym',sdftyp)) then
      call excpym(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('pyv',sdftyp)) then
      call excpyv(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('gga91',sdftyp)) then
      call excg91(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('ev',sdftyp)) then
      call excev(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      elseif(ifkey('pbe',sdftyp)) then
      call excpbe(wk(1,1),wk(1,3),wk(1,5),dr,xr,meshr,1,meshr)
      else
      call errtrp(1,'cstate','sdftyp '//sdftyp//' not serviced')
      endif
      if(ip .gt. 1 .and. abs(e-eold) .lt. tol) go to 120
c     enew=5d-1*(e+eold)
      eold=e
c     e=enew
      wk(meshr,3)=v(meshr)
      do 60 k=1,meshr-1
   60 wk(k,3)=v(k)-amagic*wk(k,3)
c     write(*,'(1x,a,i3,a,i3)')'   nclr=',nclr,'   j=',j
c     call errtrp(2,'cstate','no convergence')
  120 continue
  100 corlvl(j)=e
      call reconf(e,ebtm,config(j),rcnfg(j),tchc,tchs,fcs,ew,ez,ng,mxl
     &           ,mxlcmp,jj,sra,ls,opnc(jj),rearth)
      rorg(j)=0d0
      if(sicon) then
      do 130 k=1,meshr
  130 wk(k,5)=-wk(k,6)*wk(k,1)+amagic*wk(k,5)
      esic=esic-rcnfg(j)*fintgr(wk(1,5),dr,xr,meshr)
      endif
      anc=anc+rcnfg(j)*rin
      tint=tint+rcnfg(j)*(1d0-rin)
c     write(*,'(1x,a,i3,a,i3,3f12.5)')
c    &   ' nclr=',nclr,' j=',j,ebtm,e,rcnfg(j)
      do 140 k=1,meshr-1
  140 ro(k)=ro(k)+rcnfg(j)*wk(k,1)*xr(k)**2
      call extorg(rorg(j),wk,xr)
      if(abs(rorg(j)) .lt. 1d-20) rorg(j)=0d0
      rorg(j)=rcnfg(j)*rorg(j)
c     write(*,'(1x,a,3f12.6)')'e,esic,rin=',e,esic,rin
c     if(.not. sicon .and. abs(e-ebtm) .lt. small) then
c     write(*,'(a,i3,a,2f12.5)')
c    & '   ***msg in cstate...corelevel near ebtm found for nclr='
c    &     ,nclr,' e, ebtm=',e,ebtm
c     write(*,'(a,i3)')
c    & '   ***msg in cstate...corelevel near ebtm found for nclr=',nclr
c     ---the following part is to notify those core levels that
c        are to  be counted as valence states.
c     write(*,'(a,2i3,3f12.5)')
c    &    'j,jj,ebtm,e,config',j,jj,ebtm,e,rcnfg(j)
c     if (e .gt. ebtm .and. .not. opnc(jj) .and. jj .le. mxlcmp) then
      if (e .gt. ebtm .and. abs(rcnfg(j)) .lt. abs(config(j))) then
      config(j)=-abs(config(j))
      endif
      endif
   40 continue
c     write(*,*)
      ro(meshr)=tint
      end
