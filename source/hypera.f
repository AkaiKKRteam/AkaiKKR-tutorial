      subroutine hypera(corlu,corld,ew,ez,vu,vd,dr,xr
     &                 ,rorgu,rorgd,meshr,up,dn,sra,hhf,anclr,hforb)
c-----------------------------------------------------------------------
c     calculate the hyperfine fields. in relativistic case the brite
c     integral is calculated to obtain the hyperfine fields.
c     coded by h.akai, 1983, juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 corlu(18),corld(18),rorgu(20),rorgd(20)
     &      ,redcu(20),redcd(20),hc(18),ns(6),si(18),up(18),dn(18)
     &      ,vu(meshr),vd(meshr),dr(meshr),xr(meshr)
      character asymbl*2
      logical sra
      data ns/1,2,4,7,10,14/, bmagnt/9.274015d-21/, bohr/5.29177d-9/
      pi=4d0*atan(1d0)
      c=1d-3*(8d0*pi/3d0)*bmagnt/bohr**3
      factor=c/4d0/pi
c     --- factor=524.3044d0/4d0/pi
c     write(6,1000)
c1000 format(/'   ***msg in hypera...')
      nc=0
      do 40 j=1,18
   40 if(up(j)+dn(j) .gt. 0d0) nc=nc+1
      call relred(redcu,corlu,ew,ez,vu,dr,xr,meshr,nc,sra)
      call relred(redcd,corld,ew,ez,vd,dr,xr,meshr,nc,sra)
      shf=(rorgu(20)+rorgd(20))/4d0/pi
c
      if(nc .eq. 0) then
      hhf=factor*(rorgu(20)*redcu(16)-rorgd(20)*redcd(16))
      write(6,1100)asymbl(anclr),hhf
 1100 format(/'   hyperfine field of ',a,/3x,f10.3,' kG')
      write(6,1400)shf
 1400 format(/'   charge density at the nucleus'/3x,f14.4)
      return
      endif
c
      hfc=0d0
      sic=0d0
      do 10 j=1,nc
      rorgu(20)=rorgu(20)-rorgu(j)
   10 rorgd(20)=rorgd(20)-rorgd(j)
      hfv=factor*(rorgu(20)*redcu(16)-rorgd(20)*redcd(16))
      siv=(rorgu(20)+rorgd(20))/4d0/pi
      do 20 j=1,6
      jm=j-1
      js=ns(j)
      if(ns(j) .gt. nc) go to 30
      hc(j)=factor*(rorgu(js)*redcu(js)-rorgd(js)*redcd(js))
      si(j)=(rorgu(js)+rorgd(js))/4d0/pi
      sic=sic+si(j)
   20 hfc=hfc+hc(j)
      jm=6
   30 hhf=hfc+hfv+hforb
      write(6,1200)asymbl(anclr),hhf,hfc,hfv,hforb
      write(6,1300)(hc(j),j,j=1,jm)
      write(6,1500)shf,sic,siv
      write(6,1600)(si(j),j,j=1,jm)
 1200 format(/'   hyperfine field of ',a,/3x,f10.3,' kG'
     &       ,' (core=',f10.3,' kG','  valence=',f10.3,' kG'
     &       ,'  orbital=',f10.3,' kG )')
 1300 format
     &  ('   core contribution'
     &   /4(3x,f10.3,' kG(',i1,'s)'))
 1500 format(/'   charge density at the nucleus'/3x,f14.4
     &       ,' (core=',f14.4,'  valence=',f14.4,' )')
 1600 format
     &  ('   core contribution'
     &   /4(3x,f14.4,'(',i1,'s)'))
      return
      end
