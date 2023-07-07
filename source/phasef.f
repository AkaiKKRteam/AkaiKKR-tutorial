      subroutine phasef(v,tchc,tchs,fcs,mxl,mxlcmp,elvl,tm,ew,ez,ng
     &                ,dr,xr,meshr,rstr,wk,reltyp,is,opnc)
c-----------------------------------------------------------------------
c     --- spin-orbit and SIC correction included in this version ---
c     Construct the Tchebycheff expansion of the phase function.
c     Originally coded by H.Akai, 1983, julich
c     Revised by H.Akai, 1994, Osaka
c     Adopted to spin-orbit version, Nov 1995, Osaka
c     Modified by H. Akai, 25 Aug. 1999, Osaka
c     A separate treatment for the f resonance states is introduced
c     by H. Akai, Aug. 2014, Tokyo
c     SIC part is coorporated
c     by H. Akai 15 Jan. 2015, Tokyo
c     Adapted to open-core treatment of f-states
c     H. Akai, Tokyo, July 18, 2017
c     H. Akai, Tokyo, Jan. 24, 2018
c     Major revise
c     Add the information aoubt s^2+s^2, which is used to to calculate
c     c+is. When c+is tunrs out rather small meaning that the accuracy
c     is not enough, it is necessary to calculate c+is as
c     (c^2+s^2)/(c-is).
c     Another modification is an eraboration that make it possible
c     to handle shallow core states. This makes a smooth crossover
c     from a valence state  possible.
c     This program is now used not only f-states but also for all other
c     states.
c     by H. Akai, Tokyo, Oct.29, 2020
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v(meshr),dr(meshr),xr(meshr),wk(meshr,2*mxl+1)
     &      ,elvl(ng),tm(ng,ng),fcs(4,mxl**2),tchc(ng,mxl**2)
     &      ,tchs(ng,mxl**2),rstr(meshr,mxl**2,ng)
      real*8,allocatable::sn(:,:),cn(:,:),tc1(:),tc2(:),tc(:)
     &      ,td(:),rj(:),dj(:),rn(:),dn(:),tchr(:,:),rxg(:,:),dxg(:,:)
      complex*16 t
      character reltyp*6
      logical sra,ls,ifkey,chk,resn,opnc
      data c/274.0720442d0/, chk/.false./, jrmx/3/, itrmx/20/
     &     ,zero/1d-3/
      mxj=mxl**2
      allocate(sn(mxj,ng),cn(mxj,ng),tc1(ng),tc2(ng),tc(ng)
     &        ,td(ng),rj(mxl),dj(mxl),rn(mxl),dn(mxl),tchr(ng,mxj)
     &        ,rxg(ng,mxj),dxg(ng,mxj))
c     ls=.true. for the spin-orbit coupling included.
      ls=ifkey('ls',reltyp)
c     sra=.true. for a scalar relativistic calculation.
      sra=ifkey('sr',reltyp)
      spn=(-1d0)**(is-1)
      lmax=mxlcmp-1
      mxjcmp=mxlcmp**2
      rtin=xr(meshr-1)
      e1=ew-ez
      e2=ew+ez
      v0=v(meshr)-v(meshr-1)
      if(.not. ls) v0=0d0
      call clrarr(tchc,ng*mxj)
      call clrarr(tchs,ng*mxj)
      call clrarr(tchr,ng*mxj)
      call clrarr(wk,meshr*(2*mxl+1))
      if(chk) write(*,*)
      if(ls) then
      do 10 k=1,meshr-1
   10 wk(k,1)=v(k)*xr(k)
      call diffn(wk,wk(1,2),dr,meshr-1)
      do 20 k=1,meshr-1
   20 wk(k,2)=(wk(k,2)-v(k))/xr(k)**2
      endif
      do 30 m=-lmax,lmax
      mm=3+lmax+m
      if(ls .and. .not. sra) then
c     --- for the 'nrlls' case, the spin-orbit effects are independent
c         on energy and hence can be considered here.
      fct=spn*dble(m)/c**2
      do 40 k=1,meshr-1
   40 wk(k,mm)=v(k)+fct*wk(k,2)
      else if(.not. ls) then
c     --- either 'nrl' or 'sra' case
      do 50 k=1,meshr-1
   50 wk(k,mm)=v(k)
      endif
   30 wk(meshr,mm)=v(meshr)
      do 60 k=1,ng
c     --- construct m and e-dependent potential for which the spin-orbit
c         coupling is included ('srals' case).
c     --- nrl cases are already considered in the previous step.
      ev=elvl(k)
      if(sra) then
c     --- sra case
      ev=ev+(ev/c)**2
      if(ls) then
c     --- if both sra and ls are true, i.e. 'srals' case.
      do 70 m=-lmax,lmax
      mm=3+lmax+m
      fct=spn*dble(m)/c**2
      do 70 kk=1,meshr-1
   70 wk(kk,mm)=v(kk)+fct*wk(kk,2)/(1d0-(v(kk)-elvl(k))/c**2)**2
      endif
      endif
c
c     --- The first step is to calculate sn, cn, and tchr for
c         each elvl(k).
      call mdbsl(ev*rtin**2,lmax,rj,dj,rn,dn)
      do 60 j=1,mxjcmp
      jj=sqrt(dble(j)-5d-1)+1
      if(ls .or. j .eq. (jj-1)**2+1) then
      m=j-1-jj*(jj-1)
      mm=3+lmax+m
      fct=spn*dble(m)/c**2
      call radial(elvl(k),jj,rstr(1,j,k),rxg(k,j),dxg(k,j),v
     &           ,wk(1,mm),dr,xr,meshr,sra)
c     endif
      dxg(k,j)=dxg(k,j)+fct*v0*rxg(k,j)
      sn(j,k)=rtin**jj*(rtin*rj(jj)*dxg(k,j)-dj(jj)*rxg(k,j))
      cn(j,k)=rtin**(1-jj)*(rtin*rn(jj)*dxg(k,j)-dn(jj)*rxg(k,j))
c
      if(.false.) then
      ea=ev**(jj-1)
      t=-cn(j,k)/sn(j,k)/ea-(0d0,1d0)*csqrt(cmplx(ev,1d-3))
      write(*,'(1p,10e15.7)')elvl(k),t
      endif
c
      do 80 n=1,ng
   80 tchr(n,j)=tchr(n,j)+tm(n,k)*rxg(k,j)
      endif
   60 continue
c     j=10
c     do 22 k=1,ng
c  22 write(*,'(3f15.7)')elvl(k),cn(j,k),sn(j,k)
c  22 write(*,'(2f15.7)')elvl(k),rstr(meshr-1,j,k)
c     stop

c
c     --- Rmaining parts are performed for each (l,m).
      do 90 j=1,mxjcmp
      jj=sqrt(dble(j)-5d-1)+1
      if(ls .or. j .eq. (jj-1)**2+1) then
      m=j-1-jj*(jj-1)
      fct=spn*dble(m)/c**2
      mm=3+lmax+m
c
c     --- The second step is to fix an approximate energy that gives the
c         zero of rxg and dxg.
      le0=1
c     --- find a zero with the lowest energy.
c     --- check if rxg changes its sign between l and l-1
      do 110 l=2,ng
  110 if(rxg(l,j)*rxg(l-1,j) .lt. 0d0) le0=l
      if(le0 .gt. 1) then
c     --- Resonance for which rx=0 or dx=0 may exist.
      resn=.true.
c     --- Newton-Raphson method on a Tchebyshev expansion
c         is used to obtain e0
      e0=5d-1*(elvl(le0)+elvl(le0-1))
      itrf=0
      do 120 itr=1,itrmx
      call gntcds(e0,ew,ez,tc,td,ng)
      r1=0d0
      d1=0d0
      do 130 n=1,ng
      r1=r1+tchr(n,j)*tc(n)
  130 d1=d1+tchr(n,j)*td(n)
      de=-r1/d1
      e0=e0+de
      if(abs(de) .lt. 1d-6) exit
  120 itrf=itr
      if(itrf .eq. itrmx .or. e0 .lt. e1 .or. e0 .gt. e2) then
      call errtrp(2,'phasef','e0 not found')
      write(*,'(a,i2)') '      for j=',j
c     stop
c     --- No resonance found within the energy interval.
      resn=.false.
      er=e2
      else
      er=e0
c
c     --- determine the width of resonance.
      if(jj .gt. jrmx .and. resn) then
      x1=e0
      y1=d1
      x2=elvl(ng)
      y2=rxg(ng,j)/(x2-er)
      do 140 itr=1,10
      x3=5d-1*(x1+x2)
      call gntcs(x3,ew,ez,tc,ng)
      y3=0d0
      do 150 k=1,ng
  150 y3=y3+tchr(k,j)*tc(k)
      y3=y3/(x3-er)
      if(abs(y1+y2) .gt. 2d0*abs(y3)) then
      x2=x3
      y2=y3
      else
      x1=x3
      y1=y3
      endif
  140 continue
      delt=sqrt(3d0)*abs(x3-er)
      else
      delt=2d0*ez
      endif
      endif
      else
c     --- When no resonance exists, e0 is set to the highest
c         energy of energy interval because e0 normally
c         corresponds to the top of band.
      resn=.false.
      er=e2
      delt=2d0*ez
      endif
c
c     --- The third step is to modify cn's and sn's so as to be
c         safely fitted by Chebyshev polinomials.
      e0=er
      ev=e0
      if(sra) then
      ev=ev+(ev/c)**2
      if(ls) then
      do 160 k=1,meshr-1
  160 wk(k,mm)=v(k)+fct*wk(k,2)/(1d0-(v(k)-e0)/c**2)**2
      endif
      endif
      call mdbsl(ev*rtin**2,jj-1,rj,dj,rn,dn)
c     --- It is safer not to use values obtained through
c         Tchebyshev expansion for fcs(1,j) and fcs(2,j),
c         in particular for f-states. The derivative of
c         radial functions changes rapidly near the
c         resonance and the Tchebyshev expansion might not
c         be accurate enough to fix these biases.
      call radial(e0,jj,wk,rx,dx,v,wk(1,mm),dr,xr,meshr,sra)
      dx=dx+fct*v0*rx
c
c     --- cx and sx at e0 are used as biases.
      fcs(1,j)=rtin**(1-jj)*(rtin*rn(jj)*dx-dn(jj)*rx)
      fcs(2,j)=rtin**jj*(rtin*rj(jj)*dx-dj(jj)*rx)
      fcs(3,j)=e0
c     write(*,'(i3,f12.5)')j,e0
c
      if(jj .le. jrmx .or. .not. resn) then
c     if(.true.) then
c     --- no special treatment needed
      fcs(4,j)=0d0
      do 170 k=1,ng
      cn(j,k)=(cn(j,k)-fcs(1,j))/(elvl(k)-fcs(3,j))
  170 sn(j,k)=(sn(j,k)-fcs(2,j))/(elvl(k)-fcs(3,j))
      else
c
c     --- Suitable treatments when a resonance exists.
c     --- For details see the note 25 Mar. 2014 and 11 Apr. 2014
c         of H. Akai. The points are the following:
c         R(e), as a function of energy e, is well fitted in the
c         form of (R(e)-R(e0))/(e-e0)= A(e)*D/((e-e0)^2 + D^2),
c         where R(e) is the value of radial wave function at r=r0.
c         Now the value e0 is chosen such that R(e0)=0 is satisfied.
c         Above fitting is reasonable because for resonant states
c         R(e) sudenly changes from a nearly constant positive
c         value to a nearly constant negative value at e=e0.
c         In this case the function A(e) is expected to be a rather
c         moderate fucntion of e.
c         Similar behaviers with the same D appear also in cn and sn
c         fucntions. Thus, the fitting procedures using the same D
c         are also applied for cn and sn. In the open core
c         cases, further elaboration is needed, see below.
c
      fcs(4,j)=delt
      do 180 k=1,ng
c     write(*,'(3f12.5)') elvl(k),cn(j,k),sn(j,k)
      gam=(elvl(k)-fcs(3,j))/fcs(4,j)
      cn(j,k)=(cn(j,k)-fcs(1,j))*(gam+1d0/gam)
  180 sn(j,k)=(sn(j,k)-fcs(2,j))*(gam+1d0/gam)
      if(opnc .and. jj .eq. 4) then
c
c     --- Remove the resonant part from cn and sn
c     --- rxg and dxg are modified so as not to show
c         resonant behavior. 
c         New rxg and dxg are
c         (rxg(e)-rxg(e0)) * (e-e0)/sqrt((e-e0)**2+d**2) + rxg(e0)
c         (dxg(e)-dxg(e0)) * (e-e0)/sqrt((e-e0)**2+d**2) + dxg(e0)
c         The factor (e-e0)/sqrt((e-e0)**2+d**2) becomes sign(e-e0)
c         for |e-e0| >> d and including this factor eliminates sudden
c         jumps in rxg and dxg at the resonance, and hence, removes
c         the resonace.
c
      do 182 k=1,ng
      sgn=sqrt((elvl(k)-e0)**2+fcs(4,j)**2)/(elvl(k)-e0)
      rxg(k,j)=(rxg(k,j)-rx)*sgn+rx
      dxg(k,j)=(dxg(k,j)-dx)*sgn+dx
      ev=elvl(k)
      if(sra)ev=ev+(ev/c)**2
      call mdbsl(ev*rtin**2,lmax,rj,dj,rn,dn)
      cn(j,k)=rtin**(1-jj)*(rtin*rn(jj)*dxg(k,j)-dn(jj)*rxg(k,j))
  182 sn(j,k)=rtin**jj*(rtin*rj(jj)*dxg(k,j)-dj(jj)*rxg(k,j))
      fcs(1,j)=1d0
      fcs(2,j)=1d0
      fcs(3,j)=e2+1d0
      fcs(4,j)=0d0
      endif
c
      endif
c
c     --- The final step is to calculate the Chebyshev expansion
c         coefficints.
      do 190 l=1,ng
      do 190 k=1,ng
      tchc(l,j)=tchc(l,j)+cn(j,k)*tm(l,k)
  190 tchs(l,j)=tchs(l,j)+sn(j,k)*tm(l,k)
      endif
   90 continue
      if(.not. ls) then
      do 220 j=1,mxlcmp
      i0=(j-1)**2+1
      do 220 i=i0+1,j**2
      do 230 n=1,ng
  230 call equarr(rstr(1,i0,n),rstr(1,i,n),meshr)
      call equarr(tchc(1,i0),tchc(1,i),ng)
      call equarr(tchs(1,i0),tchs(1,i),ng)
  220 call equarr(fcs(1,i0),fcs(1,i),4)
      endif
      deallocate(sn,cn,tc1,tc2,tc,td,rj,dj,rn,dn,tchr,rxg,dxg)
      end
