      subroutine curie(wkc,ff,e,kmx,ew,ez,mxl,ng,iblk,natm,ntyp
     &               ,itype,vc,sra,a,nk,np,tch,pexf,prr,hh,ncub,clks
     &               ,last,wtkp,wt,nd,t,ncmp,ncmpx,tcpa,korder,convrg
     &               ,iatm,ls,nk3,lmxtyp,mxlcmp,msiz,lmxblk,iatmp,itblk
     &               ,its,uu,conc,length)
c-----------------------------------------------------------------------
c     This program calculates the Curie temperature of an electron
c     system by mapping it to the Heisenberg model that might be valid
c     for a small tilting of the local magnetic moments. The procedure
c     is after Liechtenstein[1].
c     [1] A.I. Liechtenstein, M.I. Katsnelson, V.P. Antropov, and
c     V.A. Gubanov, J. Magn. Magn. Mater. 67, 65 (1987)
c
c     Coded by H. Akai, 16 Nov. 2016, Tokyo
c     Copyright 2016 H. Akai. All rights reserved.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 wkc(mxl**2,mxl**2,7),ff(mxl**2,mxl**2,natm,kmx,2)
     &          ,t(mxl**2,ncmpx,kmx,2),e(kmx,2)
     &          ,tcpa(mxl**2,mxl**2,natm,kmx,2)
     &          ,tch((2*mxl-1)**2,ng,nk,nd),pexf(np,nk,nd)
     &          ,wt(ng,3,kmx,2),uu(2*mxl-1,2*mxl-1,2)
      real*8     prr(np,nk),hh(np,(2*mxl-1)**2,nk),clks(last)
     &          ,wtkp(nk),conc(ncmpx)
      integer    itype(natm),iblk(natm,natm),ncub(last)
     &          ,lmxblk(nd),itblk(5,its),ncmp(ntyp),korder(nk)
     &          ,iatm(ntyp),lmxtyp(ntyp),mxlcmp(ncmpx),iatmp(natm)
      logical    convrg(kmx),dspmf,sra,ls
      call errtrp(2,'curie','curie not implemented')
      end
