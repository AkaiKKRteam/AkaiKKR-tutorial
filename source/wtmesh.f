      subroutine wtmesh(e,wt,kmx)
c-----------------------------------------------------------------------
c     This program is used to calculate weight at
c     each mesh point if Sympson's formula is used
c     for a simple integraion along energy contour.
c     If the value of integand at e(k) is given as f(k),
c     the integral is Sum_k wt(k)*f(k).
c     coded by H. Akai, Tokyo, Dec. 3, 2021.  
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 e(kmx),wt(kmx),c1,c2,c3,e1,e2,e3
      wt(:)=(0d0,0d0)
      do 10 k=1,kmx-2,2
      e3=(e(k+2)**3-e(k)**3)/3d0
      e2=(e(k+2)**2-e(k)**2)/2d0
      e1=e(k+2)-e(k)
      c3=1d0/(e(k)-e(k+1))/(e(k)-e(k+2))
      c2=-c3*(e(k+1)+e(k+2))
      c1=c3*e(k+1)*e(k+2)
      wt(k)=wt(k)+c1*e1+c2*e2+c3*e3
      c3=1d0/(e(k+1)-e(k))/(e(k+1)-e(k+2))
      c2=-c3*(e(k)+e(k+2))
      c1=c3*e(k)*e(k+2)
      wt(k+1)=wt(k+1)+c1*e1+c2*e2+c3*e3
      c3=1d0/(e(k+2)-e(k))/(e(k+2)-e(k+1))
      c2=-c3*(e(k)+e(k+1))
      c1=c3*e(k)*e(k+1)
   10 wt(k+2)=wt(k+2)+c1*e1+c2*e2+c3*e3
      s=0d0
c     do 20 k=1,kmx
c  20 s=s+wt(k)
c     write(*,*)s
      end
