      subroutine relred(redc,corlvl,ew,ez,v,dr,xr,meshr,nc,sra)
c-----------------------------------------------------------------------
c     non-local breit integral for the hyperfine interaction strongly
c     reduce the hyperfine field obtained by the classical expression.
c     this program calculate this reduction factor by performing brite
c     integral. for valence state the central energy of the band is
c     used to calculate the factor because its energy dependence is
c     quite samll across the valence band.
c     coded by h.akai, 1983, juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 v(meshr),dr(meshr),xr(meshr),redc(20),l(15),corlvl(15)
      logical sra
      data l/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
      if(sra)then
      emx=ew-ez
      do 20 n=1,nc
      j=l(n)+1
      e=corlvl(n)
   20 call brintg(e,emx,redc(n),j,v,dr,xr,meshr)
      do 30 j=1,3
      jj=j
   30 call brintg(ew,emx,redc(15+j),jj,v,dr,xr,meshr)
      else
      do 10 i=1,20
   10 redc(i)=1d0
      endif
      end
