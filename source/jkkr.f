      subroutine jkkr(wkc,ff,e,kmx,vkp,ew,ez,mxl,ng,iblk,natm
     &               ,ntyp,itype,vc,sra,a,nk,np,tch,pexf,prr,hh,ncub
     &               ,clks,last,wtkp,wt,nd,t,tc,gfree,ess,ncmp,ncmpx
     &               ,tcpa,korder,convrg,iatm,ls,nk3,lmxtyp,mxlcmp
     &               ,msiz,lmxblk,iatmp,itblk,its,atmicp,r,g,protat
     &               ,irotat,uu,isymop,type,conc,dlimit,length)
c-----------------------------------------------------------------------
c     This program calculates J_ij parameters that map an itinerant
c     electron systems to the Heisenberg model that might be valid for
c     a small tilting of the local magnetic moments. The procedure
c     is after Liechtenstein[1].
c     [1] A.I. Liechtenstein, M.I. Katsnelson, V.P. Antropov, and
c     V.A. Gubanov, J. Magn. Magn. Mater. 67, 65 (1987)
c
c     Coded by Hisazumi Akai on 7/18/14.
c     Revised by H. Akai, 2 Dec. 2015, Tokyo
c     Copyright 2014 H. Akai. All rights reserved.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(small=1d-8,big=1d20)
      complex*16 wkc(mxl**2,mxl**2,7),ff(mxl**2,mxl**2,natm,kmx,2)
     &          ,t(mxl**2,ncmpx,kmx,2),e(kmx,2),gfree(kmx)
     &          ,ess(kmx,(2*mxl-1)**2),tcpa(mxl**2,mxl**2,natm,kmx,2)
     &          ,tch((2*mxl-1)**2,ng,nk,nd),pexf(np,nk,nd),tc(ng,kmx)
     &          ,wt(ng,3,kmx,2),uu(2*mxl-1,2*mxl-1,2)
      real*8     vkp(3,nk),prr(np,nk),hh(np,(2*mxl-1)**2,nk),clks(last)
     &          ,wtkp(nk),atmicp(3,natm),r(3,3),g(3,3)
     &          ,protat(9,48),dk(3),conc(ncmpx)
      integer    itype(natm),iblk(natm,natm),ncub(last)
     &          ,lmxblk(nd),itblk(5,its),ncmp(ntyp),korder(nk)
     &          ,iatm(ntyp),lmxtyp(ntyp),mxlcmp(ncmpx),iatmp(natm)
     &          ,irotat(natm,48),nc(3),isymop(48)
      logical    convrg(kmx),jtabl,ls,sra
      call errtrp(2,'jkkr','jkkr not implemented.')
      end
