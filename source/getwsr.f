      subroutine getwsr(wsradi,ro,dr,xr,meshr,ntyp,ncmp,ncmpx,conc
     &                 ,anclr,wk,iopt)
c-----------------------------------------------------------------------
c     Given radial charge dnesities ro's, this program gives plausible
c     values of Wigner-Seitz rasius for each type of site. There are two
c     versions:
c     1) The position where the charge neutrality inside the atomic
c        sphere is fulfilled.
c     2) The posision where ro takes a minumum in the vicinity of the
c        atomic radius defined  by the function qvolume is adopted as
c        WS radius.
c     Coded by H. Akai, Tokyo, Dec. 27, 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mpol=3,nstart=200,nh=5,ndata=2*nh+1)
      real*8 wsradi(ntyp),ro(meshr,ncmpx),dr(meshr,ncmpx)
     &      ,xr(meshr,ncmpx),conc(ncmpx),anclr(ncmpx),p(mpol),wk(meshr)
      integer ncmp(ntyp)
c
      pi=4d0*atan(1d0)
      ij=0
      do 10 i=1,ntyp
      wsradi(i)=0d0
      do 10 j=1,ncmp(i)
      ij=ij+1
      if(abs(anclr(ij)) .lt. 1d-10) then
c
c     --- in the case of a vacancy
      x=(qvolum(1d0,2)*3d0/4d0/pi)**(1d0/3d0)
c
      else if(iopt .eq. 1) then
c
c     --- method (1)
      do 20 k=1,meshr-1
   20 wk(k)=ro(k,ij)*dr(k,ij)*xr(k,ij)**2
      z=3d0*anclr(ij)
      charg=0d0
      do 30 k=1,meshr-3,2
      kk=k+2
      charg=charg+wk(k)+4d0*wk(k+1)+wk(k+2)
      if(charg .gt. z) exit
   30 continue
c     if(kk .eq. meshr-1)
c    &    call errtrp(1,'getwsr','initial wsradi too small')
      x=xr(kk,ij)-dr(kk,ij)*(charg-z)/3d0/wk(kk)
      else if(iopt .eq. 2) then
c
c     --- method (2)
c     write(*,*)'ij=',ij
c     write(*,'(2f15.5)')(xr(k,ij),ro(k,ij),k=meshr-200,meshr-1)
      do 40 k=meshr-nstart,meshr-1-nh,nh
      if(ro(k-nh,ij) .gt. ro(k,ij) .and. ro(k,ij) .lt. ro(k+nh,ij)) then
      call lsqftr(xr(k-nh,ij),ro(k-nh,ij),ndata,p,mpol)
      x=-p(2)/p(1)/2d0
      go to 50
      endif
   40 continue
      call errtrp(2,'getwsr','wsradi was not obtained')
      write(*,*)'i,j,ij=',i,j,ij
      write(*,'(a,13x,a)')'xr','ro'
      write(*,'(2f15.5)')(xr(k,ij),ro(k,ij),k=meshr-nstart,meshr-1)
      stop
   50 continue
      else
      call errtrp(1,'getwsr','illegal iop')
      endif
   10 wsradi(i)=wsradi(i)+conc(ij)*x
      end
