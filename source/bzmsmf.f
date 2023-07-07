      subroutine bzmsmf(ff,detl,tmp1,tmp2,nvec,mxl,natm,kmx,convrg
     &                 ,lmxtyp,itype)
c---------------------------------------------------------------------
c    Given tmp1 and tmp2, this routine sum up partial sums into
c    complete sums ff and detl.
c    Coded by H. Akai, Tokyo, Jan. 26, 2021
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 tmp1(mxl**2,mxl**2,natm,kmx,nvec),tmp2(kmx,nvec)
     &          ,ff(mxl**2,mxl**2,natm,kmx),detl(kmx)
      integer lmxtyp(*),itype(natm)
      logical convrg(kmx)
      do 20 ivec=1,nvec
      do 20 k=1,kmx
      if(.not. convrg(k)) then
      do 10 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      do 10 mc=1,mmx
      do 10 mr=1,mmx
c     ---add up the site-diagonal block of gi to obtain f(k,l1,l2,i).
   10 ff(mr,mc,i,k)=ff(mr,mc,i,k)+tmp1(mr,mc,i,k,ivec)
      detl(k)=detl(k)+tmp2(k,ivec)
      endif
   20 continue
      end
