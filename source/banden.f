      subroutine banden(detl,e,wt,bnde,kmx,ng)
c-----------------------------------------------------------------------
c     Calculate the band energy by integrating the phase.
c     coded by H.Akai, 1985, Juelich
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      complex*16 detl(kmx),e(kmx),wt(ng,3,kmx),a,b,c
c    &          ,d
c     p=dble(detl(kmx)-detl(1))/dble(e(kmx)-e(1))
c     q=dble(detl(1))-p*dble(e(1))
c     do 20 k=1,kmx
c  20 detl(k)=detl(k)-p*e(k)-q
      bnde=0d0
c     write(6,1100)(k,detl(k),k=1,kmx)
c1100 format((1x,i3,2f12.5))
      do 10 k=1,kmx-2,2
      a=((detl(k+1)-detl(k))/(e(k+1)-e(k))-(detl(k+2)-detl(k+1))
     &   /(e(k+2)-e(k+1)))/(e(k)-e(k+2))
      b=(detl(k+1)-detl(k))/(e(k+1)-e(k))-a*(e(k+1)+e(k))
      c=detl(k)-(a*e(k)+b)*e(k)
   10 bnde=bnde-dimag(a*(wt(1,3,k+2)-wt(1,3,k))
     &      +b*(wt(1,2,k+2)-wt(1,2,k))+c*(wt(1,1,k+2)-wt(1,1,k)))
c     bnde=bnde+dble(detl(kmx))*dimag(e(kmx))
c     --- direct integration int EdN using trapezoidal formula.
c     d=(0d0,0d0)
c     do 10 k=1,kmx-1
c     write(*,'(i3,1p4e15.6)')k,5d-1*(e(k+1)+e(k)),detl(k+1)-detl(k)
c  10 d=d+(detl(k+1)-detl(k))*(e(k+1)+e(k))
c  10 d=d+imag((detl(k+1)-detl(k)))*dble((e(k+1)+e(k)))
c  10 d=d+imag((detl(k+1)-detl(k))*(e(k+1)+e(k)))
c     bnde=5d-1*d
c     bnde=5d-1*dimag(d)
c     write(6,1000)bnde
c1000 format('   band energy for the spin state=',f10.5)
      end
