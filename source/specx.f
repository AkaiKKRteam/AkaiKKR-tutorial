                          program specx
c---------------------------------------------------------------------
c     ------------------------------------
c     --- KKR-CPA + spin-orbit version ---
c     ------------------------------------
c      KKR-CPA band structure calculation
c      originaly coded by H.Akai and M.Akai, 1980, on ACOS
c           adapted to NEC sx1&2, Aug. 1986 at Osaka
c           adapted to CRAY xmp, Aug. 1986 at Juelich
c           adapted to HITAC s820, Nov. 1989 at ISSP
c           adapted to HP 730, Julv. 1992
c                latest version on 11 Feb. 1996
c      CPA is implemented by H.Akai, 7 Sep. 1996, Osaka
c      Adapted to spin-orbit included version, 18 April, 1997, Osaka
c      XMCD calculation, 1997
c      Fixed spin moment procedure implemented, 8 Aug. 1999, Duisburg
c      Open MP, automated LMD, tilted crystal, dispersion, 2013 
c      J_ij calculation, 2016
c      Tchebychev plus Broyden, 2016
c      Tc calculation, 2017
c      Fully dynamic memory allocation by H. Akai, Tokyo, Jan. 2018.
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     parameter(ngmx=31, msr=400 ,nbuf=1000)
      parameter(ngmx=21, msr=400 ,nbuf=1000)
c     parameter(ngmx=15, msr=400 ,nbuf=1000)
c     --- nbuf should be larger than max(natm,number of components)
c
c     note: mxlmx is l_max + 1 where l_max is the maximum angular
c           momentum used in the calculation.
c           Normally, the maximum size needed for the work area
c           is determined by the maximum of 4*natmmx+ntypmx,
c           ,ncmpx*msr,5*msr, 2*mxlmx**2*ncmpmx*msex, and
c           4*nf**3/31+1. 8*nf**3 corresponds to the number of
c           k-points used within full BZ. This number, however,
c           should not be reduced to less than 7128 even if less
c           k-points in the BZ, etc. are used.
c
      complex*16,allocatable::ff(:),tmt(:),phf(:),str(:)
     &          ,ess(:),tcpa(:),phase(:),ck(:),cj(:),urotat(:),uu(:)
     &          ,e(:),wt(:),detl(:),gfree(:),tc(:)
c
      real*8     elvl(ngmx),tm(ngmx,ngmx),r(3,3),angl(3)
c
      real*8,allocatable::ancp(:),atmicp(:)
     &          ,zdmy(:),corlvl(:),dr(:),xr(:),amdlng(:),rms(:)
     &          ,clks(:),v1(:),v2(:),v3(:),ro(:),rorg(:),tchc(:),tchs(:)
     &          ,fcs(:),rstr(:),dosef(:),esic(:),sm(:),anc(:),tof(:)
     &          ,f(:),rmtd(:),q(:),cm(:),config(:),hhf(:),xmd(:)
     &          ,xanclr(:),xrmt(:),xfield(:),xconc(:)
     &          ,anclr(:),rmt(:),field(:),conc(:),hforb(:)
     &          ,wk(:)
c
      integer,allocatable::itype(:),ncub(:),iblk(:),iwtyp(:),lmxblk(:)
     &          ,match(:),irotat(:),iatm(:),mxlcmp(:)
     &          ,iatmp(:),itblk(:),lmpair(:)
     &          ,lmxtyp(:),ncmp(:),jlmxtyp(:),jncmp(:)
c
      character  go*6,file*256,brvtyp*6,reltyp*6,sdftyp*12,magtyp*6
     &          ,outtyp*6,title*2000,bzqlty*8,record*4,cwtyp*4,pmxtyp*16
     &          ,dmy*6,trmkey*16
      character,allocatable::xtype(:)*16,xatmtyp(:)*16,xatmicv(:)*24
     &          ,type(:)*16,atmtyp(:)*16,atmicv(:)*24
c
      logical    openc,ifkey
      logical,allocatable::convrg(:)
c
      data ef0/7d-1/, dex0/5d-2/, emxr/1.2d0/, msex/201/,mse0/35/
c     data ef0/7d-1/, dex0/5d-2/, emxr/1d0/, msex/201/,mse0/45/
c     data ef0/7d-1/, dex0/5d-2/, emxr/1d0/, msex/201/,mse0/65/
c     data ef0/7d-1/, dex0/5d-2/, emxr/1d0/, msex/201/,mse0/201/
     &    ,ng/ngmx/, meshr/msr/, tol/1d-6/
     &    ,ids/0/, inv/0/
c     ----------------------------------------------------
      call setomp
      call uclock(cpu)
      do 10 icmp=1,1000
      allocate(xanclr(nbuf),xrmt(nbuf),xfield(nbuf),xconc(nbuf)
     &        ,jlmxtyp(nbuf),jncmp(nbuf),xtype(nbuf)
     &        ,xatmtyp(nbuf),xatmicv(3*nbuf))
      call readin(go,file,brvtyp,a,coa,boa,alpha,beta,gamma
     &           ,edelt,ewidth,reltyp,sdftyp,magtyp,record,outtyp
     &           ,bzqlty,maxitr,pmxtyp,ntyp,xtype,jncmp,xrmt,xfield
     &           ,jlmxtyp,xanclr,xconc,natm,xatmicv,xatmtyp,r,angl
     &           ,openc,*20)
      ncmpx=0
      ncmax=0
      mxl=1
      do 30 i=1,ntyp
      mxl=max(mxl,jlmxtyp(i)+1)
      ncmax=max(ncmax,jncmp(i))
   30 ncmpx=ncmpx+jncmp(i)
      lastmx=lstcgt(mxl)
      ndmx=natm*(natm-1)+1
      lmpmx=ncmpx/2
      ids=0
      mse=ewidth*dble(mse0)/2d0
      mse=min(2*mse+1,msex)
      dmy=brvtyp
      ibrav=ibrava(trmkey('tlt',dmy))
      nf=nfqlty(bzqlty,ibrav)
      allocate(anclr(ncmpx),rmt(ntyp),field(ntyp),conc(ncmpx)
     &        ,lmxtyp(ntyp),ncmp(ntyp),type(ntyp),atmtyp(natm)
     &        ,atmicv(3*natm))
      call equarr(xrmt,rmt,ntyp)
      call equarr(xfield,field,ntyp)
      call equarr(xconc,conc,ncmpx)
      call equarx(xtype,type,ntyp)
      call equari(jncmp,ncmp,ntyp)
      call equarr(xanclr,anclr,ncmpx)
      call equari(jlmxtyp,lmxtyp,ntyp)
      call equarx(xatmtyp,atmtyp,natm)
      call equarx(xatmicv,atmicv,3*natm)
      deallocate(xanclr,xrmt,xfield,xconc,jlmxtyp,jncmp,xtype,xatmtyp
     &          ,xatmicv)
c
      if(go .eq. 'go') then
      elseif(go .eq. 'dsp') then
c      --- display the result without any iteration process.
      ids=2
      outtyp='quit'
      maxitr=1
      else if(go .eq. 'dos') then
c     --- in addition, DOS will be displayed.
      ids=1
      outtyp='quit'
      maxitr=1
c     --- also the number of the energy mesh is set to be the maximum value.
      mse=msex
c     mse=5
c     ewidth=0.1
      else if(go .eq. 'mcd' .or. go .eq. 'xmd') then
c     --- in addition, MCD will be displayed.
      ids=3
      outtyp='quit'
      maxitr=1
c     --- also the number of the energy mesh is set the maximum value.
      mse=msex
      else if(ifkey('spc',go)) then
c     --- Bloch spectral functions are output
c      --- display the result without any iteration process.
      ids=4
      outtyp='quit'
      mse=msex
      maxitr=1
      else if(go .eq. 'fsm') then
c     --- fixed spin moment procedure used
      ids=5
      else if(ifkey('j',go)) then
c     --- j_ij calculation will be done
      ids=6
      outtyp='quit'
      maxitr=1
      else if(go .eq. 'tc') then
c     --- Curie temperature calculation will be done
      outtyp='quit'
      ids=7
      maxitr=1
c     else if(go .eq. 'prm') then
c     --- suitable rtin and ew will be calculated
c     ids=8
c     mse=msex
c     outtyp='quit'
c     maxitr=1
      else if(ifkey('one',go) .or. ifkey('singl',go)
     &       .or. ifkey('sng',go)) then
c     --- one-shot calculation with Fermi level adjustment will be done
      ids=9
      maxitr=100
      else
      go='ngo'
      endif
      if(go .ne. 'ngo') then
      nwk=max(4*natm+ntyp,ncmpx*msr,6*msr,2*mxl**2*ncmpx*mse,7128
     &       ,(2*mxl+2)*msr,(2*mxl-1)**2*mxl*mse,2*(2*mxl-1)**2*mxl*24
     &       ,4*nf**3/31+1,2*mxl**4*7)
c     --- 7128=(16*441+6*1200)/2, which is used in subroutine gsdatp
      niwk=2*nwk
      allocate(ff(mxl**4*natm*mse*2),tmt(mse*mxl**2*ncmpx*2)
     &        ,phf(mse*mxl**4*ncmpx),str(mse*mxl**2*ncmpx)
     &        ,ess(mse*(2*mxl-1)**2),tcpa(mse*mxl**4*natm*2)
     &        ,phase(mse*mxl**2*ncmpx),ck(mse*mxl**4*natm)
     &        ,cj(mse*mxl**4*natm),urotat((2*mxl-1)**2*mxl*24)
     &        ,uu(2*(2*mxl-1)**2),ancp(ncmpx*2),atmicp(3*natm)
     &        ,zdmy(ncmpx),corlvl(18*ncmpx*2),dr(msr*ncmpx)
     &        ,xr(msr*ncmpx),amdlng(natm*natm),rms(ncmpx*2)
     &        ,v1(msr*ncmpx*2),v2(msr*ncmpx*2),v3(msr*ncmpx*2)
     &        ,ro(msr*ncmpx*2),rorg(20*ncmpx*2),tchc(ng*mxl**2*ncmpx)
     &        ,tchs(ng*mxl**2*ncmpx),fcs(4*mxl**2*ncmpx)
     &        ,rstr(msr*mxl**2*ng*ncmpx*2),dosef(mxl**2*ncmpx*2)
     &        ,esic(ncmpx*2),sm(ng*mxl**2*ncmpx*2),anc(ncmpx)
     &        ,tof(mxl**2*ncmpx*2),f(mxl**2*ncmpx*ng*2),rmtd(natm)
     &        ,q(ncmpx),cm(ng*mxl**2*ncmpx),config(18*ncmpx*2)
     &        ,hhf(ncmpx),xmd(ncmpx*mse*2*2),clks(lastmx)
     &        ,itype(natm),ncub(lastmx),iblk(natm*natm),iwtyp(ntyp)
     &        ,lmxblk(ndmx),match(18*ncmpx*2),irotat(natm*48),iatm(ntyp)
     &        ,mxlcmp(ncmpx),iatmp(natm),itblk(5*natm**2)
     &        ,lmpair(2*lmpmx),e(mse*2),wt(ng*3*mse*2),detl(mse*2)
     &        ,gfree(mse),tc(mse*ng),convrg(mse),hforb(ncmpx)
     &        ,wk(nwk))
c     --- generate table that returns (j,i)->jip or jip->(j,i) maps.
c         the mapping is established by calling later a subroutine
c         entry jip as 'call jip(i,j,ji)' or 'call jipinv(i,j,ji)',
      call inijip(ncmp,ntyp,ncmax)
c     --- normalize the concentration of the component atoms.
c     --- also mxlcmp is set.
      ji1=0
      ji2=0
      do 40 i=1,ntyp
      cnorm=0d0
      do 50 j=1,ncmp(i)
      ji1=ji1+1
      mxlcmp(ji1)=lmxtyp(i)+1
   50 cnorm=cnorm+conc(ji1)
      if(abs(cnorm) .lt. 1d-30) call errtrp(1,'specx','illegal conc')
      do 40 j=1,ncmp(i)
      ji2=ji2+1
   40 conc(ji2)=conc(ji2)/cnorm
      call ty2ity(atmtyp,natm,type,ntyp,itype)
      msiz=0
      do 60 i=1,natm
   60 msiz=msiz+(lmxtyp(itype(i))+1)**2
c
      call sbtime(0,0)
      call spmain
     & (wk,wk,wk,nwk,niwk,e,wt,detl,ff,tmt,phf,str,tc
     & ,anclr,rmt,ancp,field,atmicv,atmicp,elvl,tm,clks,zdmy
     & ,corlvl,dr,xr,amdlng,rms,v1,v2,v3
     & ,ro,rorg,cm,tchc,tchs,fcs,rstr,dosef,esic
     & ,sm,anc,tof,f,rmtd,q,hhf,match,itype,ncub
     & ,iblk,itblk,iwtyp,config,go,file,brvtyp,reltyp
     & ,sdftyp,magtyp,outtyp,type,atmtyp,cwtyp,title,bzqlty
     & ,record,ef0,dex0,emxr,conc,mse,ng,mxl,tol,ids,inv
     & ,a,coa,boa,alpha,beta,gamma,edelt,ewidth,maxitr,pmxtyp,xmd
     & ,ntyp,natm,meshr,ndmx,lastmx
     & ,ncmp,ncmpx,gfree,ess,tcpa,phase
     & ,convrg,urotat,uu,irotat,ck,cj,iatm,lmxtyp
     & ,mxlcmp,msiz,lmxblk,iatmp,r,lmpair,lmpmx,angl,openc,hforb)
c
      call sbtime(-1,0)
      deallocate(ff,tmt,phf,str,ess,tcpa,phase,ck,cj,urotat
     &        ,uu,ancp,atmicp,zdmy,corlvl,dr,xr,amdlng
     &        ,rms,v1,v2,v3,ro,rorg,tchc,tchs,fcs,rstr,dosef,esic,sm
     &        ,anc,tof,f,rmtd,q,cm,config,hhf,xmd,clks,itype,ncub,iblk
     &        ,iwtyp,lmxblk,match,irotat,iatm,mxlcmp,iatmp,itblk,lmpair
     &        ,e,wt,detl,gfree,tc,convrg,hforb,wk)
      call endjip
      endif
      deallocate(anclr,rmt,field,conc,lmxtyp,ncmp,type,atmtyp,atmicv)
   10 continue
   20 call uclock(cpu)
c     write(*,'(1x,a,t15,f12.2,a)')'total cpu',cpu,' sec'
      call sysinf
      call endomp
      end
