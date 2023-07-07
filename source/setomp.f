      subroutine setomp
c-----------------------------------------------------------------------
c     get the omp wall time and set number of threads
c     coded by H. Akai 10 Aug. 2015, Tokyo
c     modified by H. Akai, Tokyo, 7 June 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*16 os,host,arch
      integer*8 nstack
      data mstack/16/
c     data nsw/16/
      save tstart,nthread,mxth,os,host,arch,ncor
c     tstart=omp_get_wtime()
      nstack=mstack*10**6
      call system_clock(iclock,irate)
      tstart=dble(iclock)/dble(irate)
      ncor=numcor(os,host,arch)
      nthread=ncor
      mxth=nthread
c     if(nthread .gt. nsw) nthread=nthread-1
      call omp_set_num_threads(nthread)
c     --- in some systems, set stacksize does not reflected in the
c         actual calculation. In that case, stacksize has to be
c         specified by setting the environment variables kmp_stacksize.
c     --- if kmp_set_stacksize is available, following line can be
c         activated.
c     call kmp_set_stacksize_s(nstack)
      return
      entry thread(n)
c-----------------------------------------------------------------------
c     get the number of threads
c-----------------------------------------------------------------------
      n=nthread
      return
      entry mxthrd(n)
c-----------------------------------------------------------------------
c     get the possible maximum number of threads
c-----------------------------------------------------------------------
      n=mxth
      return
      entry rstomp(n)
c-----------------------------------------------------------------------
c     reset thread number
c-----------------------------------------------------------------------
      nthread=n
      call omp_set_num_threads(nthread)
      return
      entry sysinf
c-----------------------------------------------------------------------
c     output the system information
c-----------------------------------------------------------------------
c     ---in the case that the data are btained through c-code,
c         the null character should be transformed into a blank.
      do 10 i=1,16
      if(ichar(os(i:i)) .lt. ichar(' ')) os(i:i)=' '
      if(ichar(host(i:i)) .lt. ichar(' ')) host(i:i)=' '
   10 if(ichar(arch(i:i)) .lt. ichar(' ')) arch(i:i)=' '
      write(*,'(2a)')' OS: ',os(1:len_trim(os))
      write(*,'(2a)')' Host: ',host(1:len_trim(host))
      write(*,'(2a)')' Machine: ',arch(1:len_trim(arch))
      write(*,'(a,i0)')' numcor: ',ncor
      return
      entry endomp
c-----------------------------------------------------------------------
c     get the end wall time and output elapsed time.
c-----------------------------------------------------------------------
c     tend=omp_get_wtime()
      call system_clock(iclock,irate)
      tend=dble(iclock)/dble(irate)
      write(*,'(1x,a,t15,f12.2,a,i3,a//)')
     & 'elapsed time ',tend-tstart,' sec (',nthread,' threads)'
      end
