      program spc
c-----------------------------------------------------------------------
c     construct gnuplot batch file and submit it.
c     coded by H.Akai, 13 June 2015, Tokyo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer kcrt(25)
      character file*256,tmpfil*6,command*256,a*80,lbl(25)*16
      logical exist,typ2
c--- get input file
      file=' '
      call getarg(1,file)
      if(file .eq. ' ') then
      write(*,'(a)')'Usage: spc file'
      stop
      endif
      len=lftfil(file)
c
c     If a directory of the same name exists, inquire cannot
c     ditect existence of a file correctly. Therefore one should
c     use 'open' statement together with error return.
      open(12,file=file,status='old',iostat=ierr)
c     inquire(file=file,exist=exist)
c     if(.not. exist) then
      if(ierr .ne. 0) then
      write(*,'(a/a)')'File not found','Usage: spc file'
      stop
      endif
c---  create a temporary file
      do 10 i=0,999
      write(tmpfil,'(a,i3.3)')'tmp',i
      inquire(file=tmpfil,exist=exist)
      if(.not. exist) go to 20
   10 continue
      write(*,'(a)')'create temporary file failed'
      go to 30
   20 open(13,file=tmpfil,status='new')
      read(12,'(a)')a
      if(index(a,'header') .gt. 0) then
      typ2=.true.
      else
      typ2=.false.
      endif
      if(typ2) then
      read(12,*)a,eb,et,mse,kblst
      read(12,*)a,(kcrt(i),lbl(i),i=1,kblst)
      close(12)
      yscale=(et-eb)/dble(mse-1)
      endif
c     --- construct a gnuplot program
c     write(13,'(a)')'#!/usr/bin/gnuplot -persist'
      write(13,'(a)')'#!/usr/bin/gnuplot'
      write(13,'(2a)')'set terminal postscript landscape noenhanced ',
     &               'monochrome dashed defaultplex "Helvetica" 14'
      write(13,'(a)')'#set output "spc.eps"'
      write(13,'(a)')'set term x11 enhanced'
      write(13,'(a)')'set pm3d map'
      write(13,'(a)')'set nokey'
      write(13,'(a)')'unset tics'
      write(13,'(a)')'unset colorbox'
      write(13,'(a)')
     & 'set palette defined (0 "white", 60 "blue", 259 "blue")'
      if(typ2) then
      write(13,'(a)')'set ytics'
      write(13,'(a,i4,a)')'set xrange[0:',kcrt(kblst)-1,']'
      write(13,'(a,f8.3,a,f8.3,a)')'set yrange[',eb,':',et,']'
      write(13,'(a)')'set ylabel "Energy relative to '//
     &         'Fermi energy (Ry)" offset -2,0'
      write(13,'(a,i4,a)')'set arrow from 0,0 to ',kcrt(kblst)-1
     &              ,',0 linewidth 0 nohead front'
      do 40 i=2,kblst-1
   40 write(13,'(a,i4,a,f8.3,a,i4,a,f8.3,a)')
     &   'set arrow from ',kcrt(i)-1,',',eb,' to ',kcrt(i)-1,',',et
     &  ,' linewidth 1 nohead front'
      do 50 i=1,kblst
      if(lbl(i)(1:len_trim(lbl(i))) .eq. '(0 0 0)') then
      write(13,'(a,i4,a,f8.3,a)')
     &    'set label "{/Symbol G}" at ',kcrt(i)-1
c    &   ,',',eb-(et-eb)/50d0,' font "symbol,11" offset -0.5,-0.5 front'
     &   ,',',eb-(et-eb)/50d0,' offset -0.5,-0.5 front'
      else
      write(13,'(3a,i4,a,f8.3,a)')
     &    'set label "',lbl(i)(1:len_trim(lbl(i))),'" at ',kcrt(i)-1
     &   ,',',eb-(et-eb)/50d0
     &   ,' right rotate by 90 font "symbol,10" front'
      endif
   50 continue
      write(13,'(3a,f12.7,a,f8.3,a)')'splot "',file(1:len)
     & ,'" using ($1):($2*(',yscale,')+(',eb,')):($3)  matrix'
      else
      write(13,'(3a)')'splot "',file(1:len),' matrix'
      endif
c     write(13,'(a)')'#pause -1'
      write(13,'(a)')'pause -1'
c     call system('cat '//tmpfil)
      call system('gnuplot '//tmpfil)
   30 close(13)
      call system('rm '//tmpfil)
      end
      function lftfil(fil)
c---------------------------------------------------------------------
c     Get rid of the blank from the string. Remaining part is
c     compressed to the left. The function returns the length of
c     the non-zero part.
c     coded by H.Akai, 1986, Juelich
c     Modified by H.Akai, 2007
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character fil*(*)
      n=len(fil)
      j=0
      do 10 i=1,n
      if(fil(i:i) .eq. ' ') go to 10
      j=j+1
      fil(j:j)=fil(i:i)
   10 continue
      lftfil=j
      if(j .ge. n) return
      fil(j+1:n)=' '
      return
      end
      subroutine chleng(a,ln)
c---------------------------------------------------------------------
c     Returns the length of the non blank part of the character a.
c     coded by H.Akai
c---------------------------------------------------------------------
      character a*(*)
      n=len(a)
      do 10 i=n,1,-1
      ln=i
      if(a(i:i) .ne. ' ') return
   10 continue
      ln=0
      return
      end
