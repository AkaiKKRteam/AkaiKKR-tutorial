      program fmg
c---------------------------------------------------------------------
c     Utility program which handle the data in the file used
c     by 'cpa91' such that rearrangement of the records,
c     including the chnage in the number of the component atoms
c     and reversing of the spin direction, becomes possible.
c     Input data cards have a form like the following:
c
c     file1  1  2  3  4
c     file2  1  2 -1 -2
c
c     In this case 4 data named respectively 1, 2, 3, and 4 of
c     inout file will be output in the sequence of 1, 2, 1, 2 and
c     further the spin direction for the last two data is reversed.
c     Another example:
c
c     file1  1  2
c     file1  1  2  1  2
c
c     In this case both the input and the output files is file1, and
c     for the input only 2 data, 1 and 2, are assumed to recorded but
c     to be expanded up to 4 data in the output.
c     See inside for further detail.
c     The program uses some routines provided by 'cpasub91'.
c     coded by H. Akai, 27 aug. 1991 in Juelich
c---------------------------------------------------------------------
      parameter (lrec=1024,mxout=1000)
      implicit real*8 (a-h,o-z)
      real*8 ef1(2),rorg1(2),ef2(2),eorg2(2),dmy(50)
      real*8,allocatable:: v1(:,:,:),corlvl1(:,:,:),anclr1(:)
     &      ,v2(:,:,:),corlvl2(:,:,:),anclr2(:)
     &      ,xr1(:,:),dr1(:,:),xr2(:,:),dr2(:,:)
      integer in(mxout),ni(mxout),ip(mxout,2)
      character card*(lrec),file1*80,file2*80,key*(lrec),fmt*20
      logical exist,large
      write(fmt,'(a,i0,a)')'(a',lrec,')'
      do 10 icomp=1,1000
      read(5,fmt,end=20,err=20)card
      call getkey(card,file1)
      if(file1 .eq. ' ') stop
      call chkfil(file1,exist)
      if(.not. exist) then
      write(6,'(a,a20,a)')' file=',file1,' does not exist'
      go to 10
      endif
      large=.false.
      na=0
      do 30 i=1,mxout
      call getkey(card,key)
      if(key .ne. ' ') then
      read(key,'(i3)')in(i)
      if(abs(in(i)) .gt. 100) large=.true.
      na=na+1
      endif
   30 continue
      read(5,fmt,end=20,err=20)card
      call getkey(card,file2)
      nb=0
      do 40 i=1,mxout
      call getkey(card,key)
      if(key .ne. ' ') then
      read(key,'(i3)')ni(i)
      nb=nb+1
      endif
   40 continue
      do 50 i=1,nb
      do 60 j=1,na
      if(abs(ni(i)) .eq. abs(in(j))) then
      if(ni(i)*in(j) .lt. 0) then
      ip(i,1)=2
      ip(i,2)=1
      else
      ip(i,1)=1
      ip(i,2)=2
      endif
      ni(i)=j
      go to 50
      endif
   60 continue
      write(6,'(a)')' corresponding data not defined...skipped'
      go to 10
   50 continue
      if(large) then
      write(*,'(1x,a,2i4)')'na,nb',na,nb
      write(*,'(1x,a,(40i4))')'in   ',(in(i),i=1,na)
      write(*,'(1x,a,(40i4))')'ni   ',(ni(i),i=1,nb)
      else
      write(*,'(1x,a,2i3)')'na,nb',na,nb
      write(*,'(1x,a,(40i3))')'in   ',(in(i),i=1,na)
      write(*,'(1x,a,(40i3))')'ni   ',(ni(i),i=1,nb)
      endif
      do 70 i=1,na
   70 in(i)=0
      do 80 i=1,nb
      j=ni(i)
   80 in(j)=in(j)+(-1)**ip(i,1)
      ipara=1
      do 90 i=1,na
   90 if(in(i) .ne. 0) ipara=0
      do 100 j=1,2
  100 write(*,'(1x,a,i1,a,(40i3))')'ip',j,'  ',(ip(i,j),i=1,nb)
      write(*,'(1x)')
c
      call getfil(24,file1)
      rewind(24)
      read(24,end=20,err=20)itr,mna,msr
      if(na .gt. mna) then
      write(*,'(a,i3,a,i3)')' number of required records',na
     &     ,' is larger than number of exsisting records',mna
      stop
      endif
      allocate(v1(msr,mna,2),corlvl1(18,mna,2),anclr1(mna),xr1(msr,mna)
     &        ,dr1(msr,mna),v2(msr,mna,2),corlvl2(18,mna,2),anclr2(mna)
     &        ,xr2(msr,mna),dr2(msr,mna))
      irec=0
c     --- 1st record read in --
c     read(24)
c     write(*,*)'1st record skipped'
      rewind(24)
      read(24,end=110,err=120)itr1,ncmp24,meshr1,(anclr1(i),i=1,na)
     &         ,(((corlvl1(k,i,j),k=1,18),i=1,na),j=1,2)
     &         ,((dr1(k,i),k=1,meshr1),i=1,na)
     &         ,((xr1(k,i),k=1,meshr1),i=1,na)
     &         ,(((v1(k,i,j),k=1,meshr1),i=1,na),j=1,2)
     &         ,(ef1(i),i=1,2),er1,ew1,ez1,edelt1,dmy
      write(*,'(1x,a)')'1st record read in successfully'
      irec=1
      if(ipara .eq. 1) then
      de=5d-1*(ef1(1)-ef1(2))
      do 130 i=1,na
      v1(meshr1,i,1)=0d0
      v1(meshr1,i,2)=0d0
      do 140 k=1,18
      corlvl1(k,i,1)=corlvl1(k,i,1)-de
  140 corlvl1(k,i,2)=corlvl1(k,i,2)+de
      do 130 k=1,meshr1-1
      v1(k,i,1)=v1(k,i,1)-de
  130 v1(k,i,2)=v1(k,i,2)+de
      ef1(1)=ef1(1)-de
      ef1(2)=ef1(2)+de
      endif
c
c     --- 2nd record read in --
      read(24,end=110,err=110)itr2,ncmp24,meshr2,(anclr2(i),i=1,na)
     &         ,(((corlvl2(k,i,j),k=1,18),i=1,na),j=1,2)
     &         ,((dr2(k,i),k=1,meshr2),i=1,na)
     &         ,((xr2(k,i),k=1,meshr2),i=1,na)
     &         ,(((v2(k,i,j),k=1,meshr2),i=1,na),j=1,2)
     &         ,(ef2(i),i=1,2),er2,ew2,ez2,edelt2,dmy
      write(*,'(1x,a)')'2nd record read in successfully'
      irec=2
c     write(*,*)nb,meshr2,er2,ew2,ez2,edelt2
      if(ipara .eq. 1) then
      de=5d-1*(ef2(1)-ef2(2))
      do 150 i=1,na
      v2(meshr2,i,1)=0d0
      v2(meshr2,i,2)=0d0
      do 160 k=1,18
      corlvl2(k,i,1)=corlvl2(k,i,1)-de
  160 corlvl2(k,i,2)=corlvl2(k,i,2)+de
      do 150 k=1,meshr2-1
      v2(k,i,1)=v2(k,i,1)-de
  150 v2(k,i,2)=v2(k,i,2)+de
      ef2(1)=ef2(1)-de
      ef2(2)=ef2(2)+de
      endif
  110 continue
      close(24)
      write(*,'(1x)')
      call getfil(25,file2)
      if(irec .ge. 1) then
      write(25)itr1,nb,meshr1,(anclr1(ni(i)),i=1,nb)
     &         ,(((corlvl1(k,ni(i),ip(i,j)),k=1,18),i=1,nb),j=1,2)
     &         ,((dr1(k,ni(i)),k=1,meshr1),i=1,nb)
     &         ,((xr1(k,ni(i)),k=1,meshr1),i=1,nb)
     &         ,(((v1(k,ni(i),ip(i,j)),k=1,meshr1),i=1,nb),j=1,2)
     &         ,(ef1(i),i=1,2),er1,ew1,ez1,edelt1,dmy
      write(*,'(1x,a)')'1st record written on successfully'
c     write(*,*)nb,meshr1,er1,ew1,ez1,edelt1
      if(irec .ge. 2) then
      write(25)itr2,nb,meshr2,(anclr2(ni(i)),i=1,nb)
     &         ,(((corlvl2(k,ni(i),ip(i,j)),k=1,18),i=1,nb),j=1,2)
     &         ,((dr2(k,ni(i)),k=1,meshr2),i=1,nb)
     &         ,((xr2(k,ni(i)),k=1,meshr2),i=1,nb)
     &         ,(((v2(k,ni(i),ip(i,j)),k=1,meshr2),i=1,nb),j=1,2)
     &         ,(ef2(i),i=1,2),er2,ew2,ez2,edelt2,dmy
      write(*,'(1x,a)')'2nd record written on successfully'
c     write(*,*)nb,meshr2,er2,ew2,ez2,edelt2
      endif
      else
      write(6,'(a)')' no records available...skipped'
      endif
      close(25)
      deallocate(v1,corlvl1,anclr1,xr1,dr1,v2,corlvl2,anclr2,xr2,dr2)
   10 continue
      stop
  120 deallocate(v1,corlvl1,anclr1,xr1,dr1,v2,corlvl2,anclr2,xr2,dr2)
   20 stop
      end
      subroutine getkey(card,key)
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(nbuf=2000)
      character card*(*),ncard*(nbuf),key*(*)
      key=' '
      ne=len(card)
      if(ne .gt. nbuf) then
      write(*,'(1x,a)')'***err in getkey...buffer overflowed'
      stop
      endif
      do 10 i=1,ne
      is=i
      if(card(i:i) .ne. ' ') go to 20
   10 continue
      return
   20 ie=is
      do 30 i=is+1,ne
      if(card(i:i) .eq. ' ') go to 40
   30 ie=i
   40 key=card(is:ie)
      if(ie .lt. ne) then
      ncard=card(ie+1:ne)
      else
      ncard=' '
      endif
      card=ncard
      return
      end
      subroutine getfil(lunit,filnam)
c----------------------------------------------------------------------
c     This program opens file with unit=lunit and file=filnam.
c     When file does not exist, a new file is created and eof
c     is marked.
c     coded by H.Akai, 1988, Tokyo(ISSP)
c     revised by H.Akai, Feb. 1993, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character filnam*(*)
      logical exs,opn
      inquire(file=filnam,exist=exs,opened=opn,number=ifil)
      if(opn) close(ifil)
      open(lunit,file=filnam,form='unformatted',status='unknown')
      write(*,'(2a)') '   file to be accessed=',filnam
      rewind(24)
      if(exs) return
c      endfile(lunit)
c      rewind(lunit)
      write(*,'(a)') '   created'
      return
      end
      subroutine chkfil(file,exist)
c----------------------------------------------------------------------
c     This program check if file exists.
c     coded by H.Akai, 1988, Tokyo(ISSP)
c     revised by H.Akai, Feb. 1993, Osaka
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical exist,opened
      character file*(*)
      inquire(file=file,exist=exist,opened=opened,number=number)
      end
