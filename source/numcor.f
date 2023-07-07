      function numcor(os,host,arch)
c-----------------------------------------------------------------------
c     This function returns the number of cores of the system.
c     Coded by H. Akai, 13 Feb. 2015, Tokyo
c     Revised by H. Akai, Tokyo, Nov. 1, 2020
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character tmpfil*24,os*(*),arch*(*),host*(*)
      integer system
      logical exist
      data inicor/4/
      os='unknown'
      host='unknown'
      arch='unknown'
      numcor=inicor
c     --- create a temporary file under /tmp to avoid confliction
c         that might occur through NFS.
      call system_clock(iclock)
      write(tmpfil,'(a,i0)')'/tmp/tmp',iclock
c     write(*,*)tmpfil
      inquire(file=tmpfil,exist=exist)
      if(exist) then
      write(*,'(a)')' create temporary file failed:'
      write(*,'(a,i2,a)')' number of cores=',inicor,' is assumed.'
      return
      endif
c
      open(13,file=tmpfil)
      istat1=system('uname >'//tmpfil)
      istat2=system('uname -n >>'//tmpfil)
      istat3=system('uname -m >>'//tmpfil)
      istat4=system('getconf _NPROCESSORS_ONLN >>'//tmpfil)
      if(istat1**2+istat2**2+istat3**2+istat4**2 .ne. 0) then
      call errtrp(2,'numcor','system unknown')
      close(13)
      istat=system('rm '//tmpfil)
      if(istat .ne. 0) call errtrp(2,'numcor','rm tmpfile failed')
      return
      endif
      rewind(13)
      read(13,*)os,host,arch,numcor
      if(index(arch,'aarch64') .gt. 0) numcor=numcor-2
      close(13)
      istat=system('rm '//tmpfil)
      if(istat .ne. 0) call errtrp(2,'numcor','rm tmpfile failed')
      end
