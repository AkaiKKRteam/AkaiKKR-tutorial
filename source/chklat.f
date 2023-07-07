      subroutine chklat(ibrav,a,coa,boa,alpha,beta,gamma,vc,r,g
     &                 ,atmicv,atmicp,itype,natm,anclr,rmt,ntyp,aref
     &                 ,conc,ncmp,ncmpx,iatm,asa,fill
     &                 ,wk,iwk,corlvl,ro,dr,xr,meshr)
c----------------------------------------------------------------------
c     This program check values of 'a' and 'rmt' and if they are
c     unrealistic (either very large or small) it trys to generate
c     some more reasonable values. For 'a', experimental or MJW data
c     are used depending on a>=1d+4 or a<=1d-4, respectively.
c     coded by H.Akai, April 1992, Osaka
c     revised by H.Akai, Feb. 1993, Osaka
c     revised by H.Akai, 12 Aug, 1995, Duisburg
c     revised by H.Akai, 14 Aug, 1997, Duisburg
c     revised by H.Akai, 26 Dec, 2018, Tokyo
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 anclr(ncmpx),conc(ncmpx),rmt(ntyp),r(3,3),g(3,3)
     &      ,atmicp(3,natm),sft(3),cn(3)
     &      ,x(3),wk(*)
     &      ,corlvl(18,ncmpx),ro(meshr,ncmpx),dr(meshr,ncmpx)
     &      ,xr(meshr,ncmpx)
      real*8,allocatable::dist(:,:),wsradi(:),atvol(:)
      character atmicv(3,natm)*24
      integer itype(natm),iatm(ntyp),ncmp(ntyp)
      integer*4 iwk(*)
      integer,allocatable::num(:)
      logical asa,wsc,elbrt
      logical*1,allocatable::fixed(:)
      data zero/1d-10/, tiny/1d-3/, iop/2/,iwsrop/2/, wsc/.false./
     &    ,elbrt/.false./
      allocate(dist(ntyp,ntyp),wsradi(ntyp),atvol(ntyp),num(ntyp)
     &        ,fixed(ntyp))
      pi=4d0*atan(1d0)
      if(iop .eq. 0) call errtrp(3,'chklat','iop=0 specified')
c     --- see how often each type appears and
c     --- asign atom number for each atom type.
      do 10 i=1,ntyp
   10 num(i)=0
      do 20 i=1,natm
      iatm(itype(i))=i
   20 num(itype(i))=num(itype(i))+1
c     iatm(1)=21
c
c     --- choose type of bravais lattice if it is not given.
      if(ibrav .eq. 0) then
c     --- see which type dominates the structure.
      idom=0
      ji=0
      wmx=0
      do 30 i=1,ntyp
      do 30 j=1,ncmp(i)
      ji=ji+1
      w=conc(ji)*dble(num(i))
      if(w .gt. wmx) then
      idom=ji
      wmx=w
      endif
   30 continue
c     --- now fix the structure.
      ibrav=inqbrv(anclr(idom))
      endif
c
      call prmvec(ibrav,coa,boa,alpha,beta,gamma,vc,r,g)
      x(1)=1d0
      x(2)=boa
      x(3)=coa
c     --- a recommended value of atomic volume.
c     idata=2 for mjw values, 1 for experimental ones.
c     get average atomic volumes, AS radii, and the unit cell volume.
      idata=2
      if(a .gt. 0.99d6) idata=1
      uvol=0d0
      ji=0
      do 70 i=1,ntyp
      atvol(i)=0d0
      do 80 j=1,ncmp(i)
      ji=ji+1
   80 atvol(i)=atvol(i)+qvolum(anclr(ji),idata)*conc(ji)
      uvol=uvol+dble(num(i))*atvol(i)
      if(atvol(i) .gt. zero) then
      wsradi(i)=(atvol(i)*3d0/4d0/pi)**(1d0/3d0)
      else
      wsradi(i)=zero
      endif
   70 continue
c     --- give a lattice constant if it is not given.
      aref=(uvol/vc)**(1d0/3d0)
      if(a .lt. 1.01d-6 .or. a .gt. 0.99d6) a=aref
c
c     --- convert the data refered to a crystal axis coordinate
c     --- to those refered to the cartesian cordinate. data forms
c     --- such as  1d0, 0.5d0a, b, 3c, x, 0.333333y are allowed.
      call getcrt(atmicv,atmicp,r,x,natm)
c
c     --- Check if the primitive unit vectors actually are
c         primitive.
      call chkprm(r)
      call rgconv(r,g,vc)
c     --- Relocate all atoms inside the primitive unit cell.
      do 150 j=1,natm
      do 160 i=1,3
      cn(i)=0d0
      do 160 n=1,3
  160 cn(i)=cn(i)+atmicp(n,j)*g(n,i)
      do 150 i=1,3
      cn(i)=dble(floor(cn(i)+1d-10))
      do 150 n=1,3
  150 atmicp(n,j)=atmicp(n,j)-cn(i)*r(n,i)
c
c     --- redifine wsradi from electron density around each atom.
      if(elbrt) then
      if(wsc) then
      call wscell(r,atmicp,iatm,itype,natm,ntyp,wsradi,g)
      else
      do 170 i=1,ntyp
  170 wsradi(i)=2.5d0*wsradi(i)/a
c 170 wsradi(i)=1.3d0*wsradi(i)/a
      call gsdatp(wk,a,vc,atmicp,r,anclr,corlvl,ro,wsradi,dr,xr
     &           ,iwk,itype,natm,ntyp,meshr,ncmp,ncmpx,conc,.false.)
      call getwsr(wsradi,ro,dr,xr,meshr,ntyp,ncmp,ncmpx,conc,anclr,wk
     &           ,iwsrop)
c     --- repeat the above procedure once again to get a better
c         charge neutrality in the WS sphere.
      if(iwsrop .eq. 1) then
      do 43 i=1,ntyp
   43 wsradi(i)=wsradi(i)/a
      call gsdatp(wk,a,vc,atmicp,r,anclr,corlvl,ro,wsradi,dr,xr
     &           ,iwk,itype,natm,ntyp,meshr,ncmp,ncmpx,conc,.false.)
      endif
      endif
c     call getwsr(wsradi,ro,dr,xr,meshr,ntyp,ncmp,ncmpx,conc,anclr,wk
c    &           ,iwsrop)
c     do 45 i=1,ntyp
c  45 wsradi(i)=wsradi(i)/a
c     call getwsr(wsradi,ro,dr,xr,meshr,ntyp,ncmp,ncmpx,conc,anclr,wk
c    &           ,iwsrop)
c     ij=0
c     do 72 i=1,ntyp
c     do 72 j=1,ncmp(i)
c     ij=ij+1
c     chrgin=fintgr(ro(1,ij),dr(1,ij),xr(1,ij),meshr)
c  72 write(*,'(a,i2,2f12.5)')
c    &     'ij,anclr,total charge=',ij,anclr(ij),chrgin
c     --- grouping the type into two groups.
c         group 1: rmt is given
c         group 2: rmt is not give, i.e. rmt(i)=0
c         vol1/2 is the total volume of group 1/2 obtained from wsradi
c         vol0 is (unnormalized) volume implied by rmt
      endif
      vol0=0d0
      vol1=0d0
      vol2=0d0
      do 180 i=1,ntyp
      if(abs(rmt(i)) .gt. 1d-3) then
      vol0=vol0+dble(num(i))*(a*rmt(i))**3
      vol1=vol1+dble(num(i))*wsradi(i)**3
      else
      vol2=vol2+dble(num(i))*wsradi(i)**3
      endif
  180 continue
      vol0=vol0*4d0*pi/3d0
      vol1=vol1*4d0*pi/3d0
      vol2=vol2*4d0*pi/3d0
      red=a*(vc/(vol1+vol2))**(1d0/3d0)
      if(abs(vol0) .gt. 1d-3) then
      red0=(vol1/vol0)**(1d0/3d0)
      else
      red0=0d0
      endif
      do 190 i=1,ntyp
c     --- wsradi is redefined again such that it is
c         proportional to rmt(i) if rmt(i) is non zero;
c         else, it is proportinal to the previously defined
c         wsradi(i).
      if(abs(rmt(i)) .gt. 1d-3) then
      wsradi(i)=rmt(i)*red*red0
      else
      wsradi(i)=wsradi(i)*red
      endif
  190 continue
      if(asa .and. elbrt) then
      do 200 i=1,ntyp
  200 rmt(i)=wsradi(i)/a
      go to 210
      endif
c     write(*,'(a/(1p5e13.6))')'wsradi=',wsradi(1:ntyp)
c
c     --- muffin-tin radii are fixed in the following blocks
      do 220 i=1,ntyp
      if(rmt(i) .gt. zero) then
      fixed(i)=.true.
      else
      fixed(i)=.false.
      endif
  220 continue
      do 230 it=1,ntyp
      i=iatm(it)
c     --- sarch the shortest distance to nearby atoms for each type
c         of the atom.
      do 240 jt=1,ntyp
  240 dist(it,jt)=1d10
      do 230 ix=-1,1
      do 230 iy=-1,1
      do 230 iz=-1,1
c     --- make translation to generate atoms in different cells
      do 250 l=1,3
  250 sft(l)=dble(ix)*r(l,1)+dble(iy)*r(l,2)+dble(iz)*r(l,3)
      do 230 j=1,natm
      jt=itype(j)
c     --- we should exclude the target atom itself.
      if((.not. (ix .eq. 0 .and. iy .eq. 0 .and. iz .eq. 0
     &   .and. j .eq. i)) .and. jt .ge. it) then
      dd=0d0
      do 260 l=1,3
  260 dd=dd+(atmicp(l,j)+sft(l)-atmicp(l,i))**2
      dd=sqrt(dd)
c     if(abs(dd-3.681502E-01) .lt. 1d-4) then
c     if(abs(dd-3.254040E-01) .lt. 1d-3) then
c     if(it .eq. 2 .and. jt .eq. 5) then
c     write(*,'(1x,a,3i3)')'type,atom,type target atom=',it,j,jt
c     write(*,'(1x,a,3i3,3f10.5)')'ix,iy,iz,sft=',ix,iy,iz,sft
c     write(*,'(1x,a,1p2e13.5)')'distance,dist=',dd,dist(it,jt)
c     write(*,'(1p6e13.5)')(atmicp(l,j)+sft(l),l=1,3)
c    &                    ,(atmicp(l,i),l=1,3)
c     endif
      if(dd .lt. dist(it,jt)) then
      dist(it,jt)=dd
c     --- dist gives the minimum distance between it- and jt-th
c     --- types of atoms
      endif
      endif
  230 continue
c     write(*,'(1x,a,2i3,1p,e13.6)')
c    & (('it,jt,dist',it,jt,dist(it,jt),jt=it,ntyp),it=1,ntyp)
c     --- check consistency of the muffin-tin radii, giving warnings
c     --- and redefining them if they are not consistent.
c     --- If iop=0, do not care about such conflictions.
      if(iop .ne. 0) then
      do 270 it=1,ntyp
      do 270 jt=it,ntyp
      if(rmt(it)+rmt(jt) .gt. dist(it,jt)+zero) then
c       write(*,'(3x,a,i2,a,f10.5,a,i2,a,f10.5)')
c    &   'rmt(',it,')=',rmt(it),'  rmt(',jt,')=',rmt(jt)
        call errtrp(2,'chklat','given rmt''s conflict; reduced')
        red=dist(it,jt)/(rmt(it)+rmt(jt))
        if(iop .eq. 1) then
c     --- for iop=1, rmt's which conflict will be reduced with
c     --- keeping other rmt's unchanged.
c     --- take care when it and jt are equivalent.
          rmtit=red*rmt(it)
          rmt(jt)=red*rmt(jt)
          rmt(it)=rmtit
        else if(iop .eq. 2) then
c     --- for iop=2, all rmt's are reduced keeping the ratios
c     --- among them unchanged.
          do 280 i=1,ntyp
  280     rmt(i)=red*rmt(i)
        else
          call errtrp(1,'chklat','illegal iop specified')
        endif
      endif
  270 continue
      endif
c     --- then we try to determine the rmt's that have not yet been
c         fixed, starting from the tightest case, proceeding step by
c         step.
      do 290 itry=1,ntyp+1
      redm=1d10
      do 300 it=1,ntyp
      do 300 jt=it,ntyp
      if(.not. fixed(it) .or. .not. fixed(jt)) then
c     --- either rmt(it) or rmt(jt) is not fixed yet.
        if(fixed(it)) then
c     --- rmt(it) is already fixed but rmt(jt) is not.
          red=(dist(it,jt)-rmt(it))/wsradi(jt)
          if(red .lt. 1d-3) red=1d10
        elseif(fixed(jt)) then
c     --- rmt(jt) is already fixed but rmt(it) is not.
          red=(dist(it,jt)-rmt(jt))/wsradi(it)
          if(red .lt. 1d-3) red=1d10
        else
c     --- neither rmt(it) nor rmt(jt) is fixed.
          red=dist(it,jt)/(wsradi(it)+wsradi(jt))
        endif
        if(red .lt. redm) then
          redm=red
          itm=it
          jtm=jt
        endif
      endif
  300 continue
c     if(redm .gt. 9d9) return
      if(redm .gt. 9d9) go to 210
      if(redm .lt. tiny)
     &   call errtrp(1,'chklat','inadequate rmt specified')
      if(.not. fixed(itm)) then
      rmt(itm)=wsradi(itm)*redm
      fixed(itm)=.true.
      endif
      if(.not. fixed(jtm)) then
      rmt(jtm)=wsradi(jtm)*redm
      fixed(jtm)=.true.
      endif
  290 continue
      call errtrp(3,'chklat','fails in determining rmt')
  210 vtin=0d0
      do 310 i=1,ntyp
  310 vtin=vtin+dble(num(i))*rmt(i)**3
      fill=(vtin*4d0*pi/3d0)/vc
      if(asa) then
      if(iop .ne. 0) then
      red=fill**(-1d0/3d0)
      fill=1d0
c     write(*,*)'vc=',vc
c     write(*,*)'sum of vtin=',vtin*4d0*pi/3d0
c     write(*,*)'red=',red
      do 320 i=1,ntyp
  320 rmt(i)=rmt(i)*red
      endif
      endif
c     do 330 i=1,ntyp
c 330 write(*,'(3x,a,i2,a,2f12.7)')'rmt(',i,')=',rmt(i),rmt(i)/wsradi(i)
c     write(*,'(3x,a,f12.7)')'volume filling=',fill
      deallocate(dist,wsradi,atvol,num,fixed)
c     stop
      end
