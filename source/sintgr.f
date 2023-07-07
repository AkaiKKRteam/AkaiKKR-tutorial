      function sintgr(f,meshr)
c-----------------------------------------------------------------------
c     Return Simpson type integral of f.
c     coded by H.Akai
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 f(meshr)
      s1=(f(1)+f(meshr-1))/2d0
      do 10 k=3,meshr-3,2
   10 s1=s1+f(k)
      s2=0d0
      do 20 k=2,meshr-2,2
   20 s2=s2+f(k)
      sintgr=(2d0*s1+4d0*s2)/3d0
      end
