      function lstcgt(mxl)
c-----------------------------------------------------------------------
c     Obtaine necessary value for 'lstcgt'
c     ---------------------------------------------------------
c      mxl    1  2  3   4    5    6    7     8     9     10
c     ---------------------------------------------------------
c      lstcgt 2 37 243 964 2854 6901 14723 28462 50848 85535
c     ---------------------------------------------------------
c     Based on the code by M.Akai, 1979, Osaka.
c     Coded by H. Akai based on the code by M.Akai (1979, Osaka),
c     Tokyo, 2018 April.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer lstmxd(10)
      data lstmxd/2,37,243,964,2854,6901,14723,28462,50848,85535/
      lstcgt=0
      if(mxl .ge. 1 .and. mxl .le. 10) then
      lstcgt=lstmxd(mxl)
      return
      elseif(mxl .gt. 10) then
      lmax=mxl-1
      ll=(lmax+1)**2
      lj=(2*lmax+1)**2
      lstcgt=0
      do 20 mr=1,ll
      call subscr(l1,m1,mr)
      do 20 mc=1,ll
      call subscr(l2,m2,mc)
      do 10 mj=1,lj
      mjj=mj
      call subscr(l3,m3,mjj)
      c=cg(l1,l2,l3,m1,m2,m3)
   10 if(abs(c) .gt. 1d-7) lstcgt=lstcgt+1
   20 lstcgt=lstcgt+1
      return
      else
      call errtrp(1,'lstcgt','illegal mxl')
      endif
      end
