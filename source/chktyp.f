      subroutine chktyp(irotat,itype,natm)
c-----------------------------------------------------------------------
c     Given the table 'irotat', which describes the rotation of atoms
c     in the unit cell, and 'itype', which gives the type of each atom,
c     this routine tells if the information given by 'itype' is
c     compatible with 'irotat'. This is useful when the unit cell
c     contains many atom. It checks whether the type assigned to each
c     atom is consistent with the crystal symmetry.
c     Coded by H. Akai, 1 Sept. 1999, Osaka.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer irotat(natm,48),itype(natm)
      logical*1,allocatable::ini(:)
      allocate(ini(natm))
      do 10 i=1,natm
   10 ini(i)=.true.
      do 20 i=1,natm
      if(ini(i)) then
      it=itype(i)
      do 30 j=i+1,natm
      if(itype(j) .eq. it) then
      do 40 l=1,48
      if(irotat(i,l) .eq. j) then
      ini(j)=.false.
      go to 30
      endif
   40 continue
      call errtrp(2,'chktyp','check atomic position or atom type.')
      write(*,'(3x,a,i3,a,i3,a/)')'atoms',i,' and',j,
     &  ' are assigned to the same type, but it obviously is wrong'
      endif
   30 continue
      endif
   20 continue
      deallocate(ini)
      end
