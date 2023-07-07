      subroutine bzmrtf(f,g,msiz,natm,irotat,urotat,isymop,kmx
     &           ,convrg,mxl,snor,lmxtyp,itype,iatmp,ls)
c-----------------------------------------------------------------------
c     Rotate the full  g(msiz,msiz) according to
c     the symmetry table given by urotat, irotat and isymop, and add up
c     the resulting matrix to g.
c     The program applicable only when rotated position staies inside
c     the orpginal unit cell, or the Bloch phase facotr exp(iRk) does
c     not apprar in g such as the case where g=|f|**2.
c
c     urotat: rotation matrices for any 24 cubic or 12 hexgonal
c             operations (inversion is omitted).
c     irotat: associated rotation in the atomic position. If atom 2
c             comes to the position where atom 1 was sitting before the
c             rotation, then irotat(1)=2, etc.
c     isymop: if isymop(iop)=1, then rotation specified by iop is
c             allowed from the symmetry.
c
c     Coded by H. Akai, 21 Nov. 2016
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(ismx=1)
c     --- ismx=1 except in the cases of Dirac or non-colinear magnetism.
      complex*16 g(ismx*msiz,ismx*msiz,kmx)
     &          ,u,urotat((2*mxl-1)**2,mxl,24)
     &          ,f(ismx*mxl**2,ismx*mxl**2,natm,kmx)
      complex*16,allocatable::wk1(:,:),wk2(:,:)
      integer irotat(natm,24,2),isymop(24,2),itype(natm),lmxtyp(*)
     &       ,iatmp(natm)
      logical convrg(kmx)
      data zero/1d-10/
      allocate(wk1(msiz,msiz),wk2(msiz,msiz))
      write(*,'(1x,24i2)')irotat
      do 10 is1=0,ismx-1
      do 10 is2=0,ismx-1
      do 10 k=1,kmx
c     --- calculation is needed only when we have not yet gotten
c         a converged results.
      if(.not. convrg(k)) then
      call clrarc(wk1,msiz**2)
      weight=0d0
      do 20 ip=1,2
      do 20 iop=1,24
c     --- only rotations compatible with the crystal symmetry
c         are considered.
      if(isymop(iop,ip) .eq. 1) then
      weight=weight+1d0
      call clrarc(wk2,msiz**2)
c     --- we first calculate wk2=f*U^+, where U^+ is the Hermite
c         conjugate of the rotation matrix U.
      do 30 l=1,mxl
      p=(-1d0)**((ip-1)*(l-1))
      mx=2*l-1
      lb=(l-1)**2
      do 30 m1=1,mx
      m0=mx*(m1-1)
      do 30 m2=1,mx
      m2m1=m0+m2
      if(abs(urotat(m2m1,l,iop)) .gt. zero) then
c     --- only non-zero elements of the rotation matrix are
c         taken into account. In order to do this procedure
c         efficiently, the order of summation appearing in
c         the matrix product is changed from the straight forward
c         way.
c     --- take account of the effects of the inversion.
      u=conjg(p*urotat(m2m1,l,iop))
      mr=lb+m1
      mc=lb+m2
      do 40 j=1,natm
      if(l .le. lmxtyp(itype(j))+1) then
      mcj=iatmp(j)+mc-1
      mrj=iatmp(j)+mr-1
      do 50 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      do 50 m=1,mmx
      mi=iatmp(i)+m-1
   50 wk2(mi,mcj)=wk2(mi,mcj)+g(mi+msiz*is1,mrj+msiz*is2,k)*u
      endif
   40 continue
      endif
   30 continue
c     --- next, we calculate U*wk2 and the results are accumulated
c         on wk1.
      do 60 l=1,mxl
      p=(-1d0)**((ip-1)*(l-1))
      mx=2*l-1
      lb=(l-1)**2
      do 60 m1=1,mx
      m0=mx*(m1-1)
      do 60 m2=1,mx
      m2m1=m0+m2
      if(abs(urotat(m2m1,l,iop)) .gt. zero) then
      u=p*urotat(m2m1,l,iop)
      mc=lb+m1
      mr=lb+m2
      do 70 i=1,natm
      if(l .le. lmxtyp(itype(i))+1) then
c     --- we also have to take it account that the atom positions
c         are rotated.
      ir=irotat(i,iop,ip)
      mri=iatmp(i)+mr-1
      mci=iatmp(ir)+mc-1
      do 80 j=1,natm
      jr=irotat(j,iop,ip)
      mmx=(lmxtyp(itype(j))+1)**2
      do 80 m=1,mmx
      mj=iatmp(j)+m-1
      mrj=iatmp(jr)+m-1
   80 wk1(mri,mj)=wk1(mri,mj)+u*wk2(mci,mrj)
      endif
   70 continue
      endif
   60 continue
      if(k .eq. 1) write(*,'(1x,1p2e14.6)')wk2(1,19)
      endif
   20 continue
c     --- now the allowed rotations have been taken, and we
c         can replace the old f by a new one. Also the summation
c         with respect to -k, which was omitted from the BZ mesh
c         for the cases of
c           i) time reversal symmetry exists
c           ii) inversion symmetry exists
c         are now taken into account.
c         "snor" is a presumed normalization factor.
      if(ls .eq. 0) then
c     --- case i) time resersal symmetry exists
      weight=2d0*weight
      call equarc(wk1,wk2,msiz**2)
      do 90 mc=1,msiz
      do 90 mr=1,msiz
   90 wk1(mr,mc)=wk2(mr,mc)+wk2(mc,mr)
      else if(isymop(1,2) .eq. 1) then
c     --- case ii) inversion symmetry exists
      weight=2d0*weight
      call equarc(wk1,wk2,msiz**2)
      do 100 i=1,natm
      mmxi=(lmxtyp(itype(i))+1)**2
      ir=irotat(i,1,2)
      do 100 im=1,mmxi
      l1=lindx(im)-1
      mr=iatmp(i)+im-1
      mrr=iatmp(ir)+im-1
      do 100 j=1,natm
      mmxj=(lmxtyp(itype(j))+1)**2
      jr=irotat(j,1,2)
      do 100 jm=1,mmxj
      l2=lindx(jm)-1
      mc=iatmp(j)+jm-1
      mcr=iatmp(jr)+jm-1
  100 wk1(mr,mc)=wk2(mr,mc)+wk2(mrr,mcr)*(-1d0)**(l1+l2)
      endif
c     --- do nothing if neither case i) nor case ii)
      factor=1d0/weight/snor
      do 110 mc=1,msiz
      do 110 mr=1,msiz
  110 g(mr+msiz*is1,mc+msiz*is2,k)=factor*wk1(mr,mc)
      do 120 i=1,natm
      mmx=(lmxtyp(itype(i))+1)**2
      mb=iatmp(i)-1
      do 120 mc=1,mmx
      do 120 mr=1,mmx
  120 f(mr,mc,i,k)=g(mb+mr,mb+mc,k)
      endif
   10 continue
c     ia=iatmp(1)-1
c     ib=iatmp(2)-1
c     write(*,*)
c     write(*,1000)((dble(g(mr,mc,1)),mc=ia+1,ia+9),mr=ib+1,ib+9)
c     write(*,*)
c     write(*,1000)((dimag(g(mr,mc,1)),mc=ia+1,ia+9),mr=ib+1,ib+9)
c1000 format(1x,1p9e14.6)
      deallocate(wk1,wk2)
      end
