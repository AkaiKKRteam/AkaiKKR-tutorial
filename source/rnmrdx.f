      function rnmrdx(b,m,ierr)
c-----------------------------------------------------------------------
c     Given character string b of length m, this program returns
c     numerical data.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character b*(m)
      read(b,*,iostat=ierr)rnmrdx
      end
