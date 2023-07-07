      subroutine pltmrk(convrg,mark,ne)
c-----------------------------------------------------------------------
c     Given logocal data 'convrg' of size ne, this program plots '*'
c     mark at the k-th column if convrg(k) is not true.
c     Coded by H. Akai, Sep. 30, 2021
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character mark*(ne)
      logical convrg(ne),notyet
      mark=' '
      notyet=.false.
      do 10 k=1,ne
      if(.not. convrg(k)) then
      mark(k:k)='*'
      notyet=.true.
      endif
   10 continue
      if(notyet) write(*,'(1x,a)')mark
      end
