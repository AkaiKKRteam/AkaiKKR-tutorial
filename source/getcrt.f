      subroutine getcrt(atmicv,atmicp,r,x,natm)
c-----------------------------------------------------------------------
c     Given character data atmicv(3,ndata) indicating atomic positions,
c     this program returns the corresponding cartesian coordinate
c     atmicp(3,natm).
c     Coded by H. Akai, Tokyo, Jan. 27, 2020
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 atmicp(3,natm),r(3,3),x(3),cartes(3)
      character atmicv(3,natm)*24,buff*80
c
c     --- convert the data refered to a crystal axis cordinate
c     --- to those refered to the cartesian cordinate. data forms
c     --- such as  1d0, 0.5d0a, b, 3c, x, 0.333333y are allowed.
      ia=ichar('a')
      iz=ichar('z')
      do 90 j=1,natm
      do 100 i=1,3
  100 cartes(i)=0d0
      do 110 i=1,3
      buff=atmicv(i,j)
c     ---get the position of the last non-blank character
      call chleng(buff,l)
      ibuff=ichar(buff(l:l))
c     ---check if the last character is a, b, c, x, y, or z.
      if(ibuff .ge. ia .and. ibuff .le. iz) then
c     ---the data is a multiple of the primitive vectors a, b, c
c     ---or the vectors x=(1,0,0), (b/a)y=(0,b/a,0), and
c        (c/a)z=(0,0,c/a).
      buff(l:l)=' '
      if(l .eq. 1) then
c     ---containing no numeric data. unity is assumed.
      atmicp(i,j)=1d0
      else
      atmicp(i,j)=redata(buff)
      endif
c
      if(ibuff .le. ia+2) then
c     ---it is a, b, or c.
      ibuff=ibuff-ia+1
      do 120 l=1,3
  120 cartes(l)=cartes(l)+atmicp(i,j)*r(l,ibuff)
      else if(ibuff .ge. iz-2) then
c     ---it is x, y, or z.
      ibuff=ibuff-iz+3
      cartes(ibuff)=cartes(ibuff)+atmicp(i,j)*x(ibuff)
      else
c     ---it is neither of them.
      call errtrp(0,'chklat','illegal character appears')
      write(*,'(1x,a)') buff(l:l)
      stop
      endif
c
      else
c     ---the data is a number, meaning multiple of (1,0,0), etc.
      atmicp(i,j)=redata(buff)
      cartes(i)=cartes(i)+atmicp(i,j)
      endif
  110 continue
      do 90 i=1,3
   90 atmicp(i,j)=cartes(i)
      end
