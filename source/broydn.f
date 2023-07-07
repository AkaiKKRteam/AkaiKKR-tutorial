      subroutine broydn(v1,v2,fm,dr,xr,u,p,g,alpha,rms,cnvq
     &                 ,meshr,ncmpx,itr)
c-----------------------------------------------------------------------
c     Implemented here according to notes of S.B.
c     Broyden iteration scheme following the papers of
c     Srivastava, J. Phys., 17, 1317 (1984),
c     R.G. Broyden, Math. Comp., 19, 577 (1965),
c     R.G. Broyden, ibid, 21, 368 (1967).
c     The method has been generalized to include a metric. The
c     definition of the necessary inner products are similar to the
c     discription given in the notes of M.Weinert. The algorithm
c     discribed in the paper Srivastava has been simplified
c     ( see notes of S.B.).
c     The files ui,vi are stored on high speed ssd memory.
c     Broyden's update treats charge and spin on the same footing
c     S. Bluegel, KFA, May 1987
c
c     Adapted to kkr-cpa program
c     by H. Akai, KFA, June 1987
c     Modified according to notes by H.A.
c     by H. Akai, Osaka, April 1988
c     Minor revision
c     by H. Akai, Tokyo, 20 Feb. 2015
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v1(meshr,ncmpx,2),v2(meshr,ncmpx,2),xr(meshr,ncmpx)
     &      ,dr(meshr,ncmpx),rms(ncmpx,2),fm(meshr,ncmpx,2)
     &      ,g(meshr,ncmpx),u(meshr,ncmpx,2,2),p(meshr,ncmpx,2,2)
      save cdiag,sdiag
c     data amix/0.012d0/
      data amix/0.05d0/
c     data amix/0.1d0/
c
      cnvq=0d0
      do 10 is=1,2
      do 10 i=1,ncmpx
      call rmserr(v1(1,i,is),v2(1,i,is),rms(i,is),dr(1,i),xr(1,i),meshr)
      if(rms(i,is) .lt. 1d-30) rms(i,is)=1d-30
   10 cnvq=max(cnvq,rms(i,is))
      call rrcnv(v1,meshr*ncmpx,meshr*ncmpx)
      call rrcnv(v2,meshr*ncmpx,meshr*ncmpx)
      do 20 is=1,2
      do 20 i=1,ncmpx
      do 20 k=1,meshr
      v1(k,i,is)=5d-1*v1(k,i,is)
   20 v2(k,i,is)=5d-1*v2(k,i,is)
c
c---->  coming block is activated only one iteration before
c       broyden iteration scheme is used
c---->  set up of : fm=fm(1)=f(rho(1))-rho(1) ;
c       metric  g := r*r*drdi (=radial volume element)
c       and generate new input potential v1
c
      if(itr .eq. 1) then
c
      do 30 is=1,2
      do 30 i=1,ncmpx
      do 30 k=1,meshr
      fm(k,i,is)=v2(k,i,is)-v1(k,i,is)
   30 v2(k,i,is)=v1(k,i,is)+amix*fm(k,i,is)
      call rrcnv(v1,meshr*ncmpx,meshr*ncmpx)
      call rrcnv(v2,meshr*ncmpx,meshr*ncmpx)
      cdiag=-alpha
c     sdiag=-8d0*alpha
      sdiag=-6d0*alpha
      do 40 i=1,ncmpx
c
c----> metric for constant shift
c
      g(meshr,i)=(xr(meshr,i)**3-xr(meshr-1,i)**3)/3d0
c
c----> metric
c
      do 40 k=1,meshr
   40 g(k,i)=xr(k,i)**2*dr(k,i)
c
      return
      endif
c
c----> v2(m)= f(rho(m))-rho(m) of all mt-spheres
c
      do 50 is=1,2
      do 50 i=1,ncmpx
      do 50 k=1,meshr
      v2(k,i,is)=v2(k,i,is)-v1(k,i,is)
c
c----> calculate dfm=f(m)-f(m-1)
c
   50 fm(k,i,is)=v2(k,i,is)-fm(k,i,is)
c
c----> loop to generate u(m)=u(k,itr)
c
      do 60 i=1,ncmpx
      do 60 k=1,meshr
      u(k,i,1,2)=-cdiag*v2(k,i,1)
   60 u(k,i,2,2)=-sdiag*v2(k,i,2)
      rewind(99)
c
      if(itr .ge. 3) then
      do 70 it=2,itr-1
      read(99)(((u(k,i,is,1),k=1,meshr),i=1,ncmpx),is=1,2)
     &       ,(((p(k,i,is,1),k=1,meshr),i=1,ncmpx),is=1,2)
      am=0d0
      do 80 is=1,2
      do 80 i=1,ncmpx
      do 80 k=1,meshr
   80 am=am+p(k,i,is,1)*v2(k,i,is)
      do 70 is=1,2
      do 70 i=1,ncmpx
      do 70 k=1,meshr
   70 u(k,i,is,2)=u(k,i,is,2)-am*u(k,i,is,1)
      endif
c
c-------->     Broyden's second method
c
c----> calculate p(m) ; convoluted with the metric g
c
      vnorm=0d0
      do 90 is=1,2
      do 90 i=1,ncmpx
      do 90 k=1,meshr
      p(k,i,is,2)=g(k,i)*fm(k,i,is)
c
c----> calculate vnorm and normalize p(m)
c
   90 vnorm=vnorm+p(k,i,is,2)*fm(k,i,is)
      vnorm=1d0/vnorm
      do 100 is=1,2
      do 100 i=1,ncmpx
      do 100 k=1,meshr
      p(k,i,is,2)=vnorm*p(k,i,is,2)
c
c----> update f(m-1)=f(m)  ; v(m)=v(m-1)
c
  100 fm(k,i,is)=v2(k,i,is)
c
c----> write on ssd u(k,itr) and p(k,itr)
c
      write(99)(((u(k,i,is,2),k=1,meshr),i=1,ncmpx),is=1,2)
     &        ,(((p(k,i,is,2),k=1,meshr),i=1,ncmpx),is=1,2)
c
c----> calculate cmm
c
      cmm=1d0
      do 110 is=1,2
      do 110 i=1,ncmpx
      do 110 k=1,meshr
  110 cmm=cmm-p(k,i,is,2)*fm(k,i,is)
c
c----> update v(m+1)
c
      do 120 is=1,2
      do 120 i=1,ncmpx
      do 120 k=1,meshr
  120 v2(k,i,is)=v1(k,i,is)+cmm*u(k,i,is,2)
      call rrcnv(v1(1,1,1),meshr*ncmpx,meshr*ncmpx)
      call rrcnv(v2(1,1,1),meshr*ncmpx,meshr*ncmpx)
      end
