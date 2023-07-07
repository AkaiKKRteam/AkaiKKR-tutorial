      function gintgr(g,dr,meshr)
c-----------------------------------------------------------------------
c     Return radial integral of g. In this program r**2 is supposed
c     to be already included in g.
c     coded by H.Akai
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 g(meshr),dr(meshr)
      s1=(g(1)*dr(1)
     &       +g(meshr-1)*dr(meshr-1))/2d0
      do 10 k=3,meshr-3,2
   10 s1=s1+g(k)*dr(k)
      s2=0d0
      do 20 k=2,meshr-2,2
   20 s2=s2+g(k)*dr(k)
      gintgr=(2d0*s1+4d0*s2)/3d0
      end
