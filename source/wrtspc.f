      subroutine wrtspc(spctrl,vkp,is,e,mse,vcrt,kcrt,kblst
     &                 ,nk1,nk3,ef,unit,iwrt)
c-----------------------------------------------------------------------
c     Given array data {{spctrl(k,kp), k=1,mse},kp=1,nk3} and
c     {e(k),k=1,mse}, this program write down plot data for the
c     dispersion relation E(k). The array {kcrt(k),k=1,kblst}
c     contains the pointer for the k-points that are on a symmetry
c     point. The data are output on a file unit=27 and 28 for spin
c     up and down, respectively.
c
c     Coded by H. Akai, 6 June 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 spctrl(mse,nk3),e(mse,2)
      real*8 vkp(3,nk1+nk3),vcrt(3,kblst),ef(2)
      integer kcrt(kblst)
      character lbl(25)*16,a*4,s(3)
      data small/1d-5/,s/'x','h','z'/
      if(iwrt .eq. 1) then
c     --- output with format (a)
      sftfct=5d-1*unit
      estep=dble(e(2,is))-dble(e(1,is))
c     write(*,'(1x,a,i2)')'#  A(E,k) for spin =',is
      write(26+is,'(a,i2)')'#  A(E,k) for spin =',is
      do 10 kp=1,nk3
      kk=nk1+kp
c     write(*,'(/1x,a,3f12.5)')'   k=',(vkp(i,kk),i=1,3)
c     write(27,'(/a,3f12.5)')'#  k=',(vkp(i,kk),i=1,3)
      if(kp .eq. 1) then
      dist=0d0
      else
      dist=dist+sqrt( (vkp(1,kk)-vkp(1,kk-1))**2
     &          +(vkp(2,kk)-vkp(2,kk-1))**2
     &          +(vkp(3,kk)-vkp(3,kk-1))**2 )*sftfct
      endif
c--- output in (x,y,z) form ----
      write(26+is,*)
      do 10 k=2,mse
   10 write(26+is,'(1x,3f15.7)')dist,5d-1*dble(e(k,is)+e(k-1,is))-ef(is)
     &    ,max(0d0,dimag((spctrl(k,kp)
     &    -spctrl(k-1,kp))/(e(k,is)-e(k-1,is))))
      return
c
      else if(iwrt .eq. 2) then
c     --- output with format (b)
c--- output in matrix form (for Igor)----
c    horizontal: energy, vertical: k-point
      do 20 kp=1,nk3
   20 write(26+is,'(1x,500f14.7)')(max(0d0,dimag((spctrl(k,kp)
     &    -spctrl(k-1,kp))/(e(k,is)-e(k-1,is)))),k=2,mse)
      return
c
      else if(iwrt .eq. 3) then
c     --- output with format (c)
c--- output in matrix form (for gnuplot)----
c    horizontal: k-point, vertical: energy
c     nkd=max(3,nk3/60)
      estep=dble(e(2,is))-dble(e(1,is))
      eb=dble(e(1,is))-ef(is)+5d-1*estep
      et=dble(e(mse,is))-ef(is)-5d-1*estep
      do 50 n=1,kblst
      kk=nk1+kcrt(n)
      lbl(n)=''
      do 50 i=1,3
      a=s(i)
c     vk=vkp(i,kk)*rta(i)
      vk=vcrt(i,n)
      if(abs(vk) .le. small)a='0'
      if(abs(vk-1d0) .le. small)a='1'
      if(abs(vk+1d0) .le. small)a='-1'
      if(abs(vk-5d-1) .le. small)a='1/2'
      if(abs(vk+5d-1) .le. small)a='-1/2'
      if(abs(vk-1.5d0) .le. small)a='3/2'
      if(abs(vk+1.5d0) .le. small)a='-3/2'
      if(abs(vk-1d0/3d0) .le. small)a='1/3'
      if(abs(vk+1d0/3d0) .le. small)a='-1/3'
      if(abs(vk-3d0/4d0) .le. small)a='3/4'
      if(abs(vk+3d0/4d0) .le. small)a='-3/4'
   50 lbl(n)=lbl(n)(1:len_trim(lbl(n)))//' '//a
      do 60 n=1,kblst
   60 lbl(n)='''('//lbl(n)(2:len_trim(lbl(n)))//')'''
c  60 if(lbl(n) .eq. '''(0 0 0)''') lbl(n)='''Gamma'''
      write(26+is,'(a/a,2e15.7,2i4/a,25(i4,1x,a))')
     &      '### header for format (c)'
     &     ,'#',eb,et,mse-1,kblst
     &     ,'#',(kcrt(i),lbl(i),i=1,kblst)
      write(26+is,'(a)')'### end of header'
      do 30 k=1,mse-1
      do 40 kp=1,nk3
   40 spctrl(k,kp)=max(0d0,dimag((spctrl(k+1,kp)
     &              -spctrl(k,kp))/(e(k+1,is)-e(k,is))))
   30 write(26+is,'(1x,500f14.7)')(dble(spctrl(k,kp)),kp=1,nk3)
      return
      else
      call errtrp(1,'wrtspc','illegal iwrt')
      endif
      end
