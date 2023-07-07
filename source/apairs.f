      subroutine apairs(d,atmicp,r,g,dlimit,protat,iatm,itype,ncmp
     &                  ,ncmpx,ntyp,natm,irotat,ipair,npairx,npair
     &                  ,ncelx,ncel,nc,jtbl,mjtblx,mjtbl,icel,isymop)
c-----------------------------------------------------------------------
c     This program genrates independent pairs of atoms, one of which
c     resides inside the unit cell centered on the origin, the other
c     being in either the same cell or a different unit cell. The
c     information created in this program includes component number,
c     site number, and the position vectors of the cell.
c
c     Coded by Hisazumi Akai on 7/18/14.
c     Revised by H. Akai, 2 Dec. 2015, Tokyo
c     Copyright 2014 H. Akai. All rights reserved.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 d(3,ncelx),atmicp(3,natm),r(3,3),g(3,3),protat(9,48),cn(3)
     &      ,s(3,2),t(3,2),s2(3)
      real*8, allocatable:: dval(:),a(:,:),aa(:,:)
      integer iatm(ntyp),irotat(natm,48),ipair(4,npairx)
     &       ,icn(3),ncmp(ntyp),itype(natm),nc(3),iswp(6),itgt(2)
     &       ,jtbl(4,npairx),icel(0:3,ncelx),icl(3),isymop(48)
      integer,allocatable:: lookup(:,:,:,:),map(:),imap(:),invrot(:,:)
      data small/1d-10/
c
c     --- memory allocation
      allocate(lookup(-nc(1):nc(1),-nc(2):nc(2),-nc(3):nc(3),natm)
     &        ,dval(npairx),map(npairx),imap(npairx),a(3,ncelx)
     &        ,aa(3,ncelx),invrot(natm,48),stat=ierr)
      if(ierr .ne. 0) call errtrp(1,'apairs','allocation fails')
c
c     --- Make up a table of neighboring cells.
      ncl=0
      do 10 n1=0,2*nc(1)
      i1=mod(n1+nc(1),2*nc(1)+1)-nc(1)
      do 10 n2=0,2*nc(2)
      i2=mod(n2+nc(2),2*nc(2)+1)-nc(2)
      do 10 n3=0,2*nc(3)
      i3=mod(n3+nc(3),2*nc(3)+1)-nc(3)
      ncl=ncl+1
      if(ncl .gt. ncelx) call errtrp(1,'apairs','ncelx too small')
c     --- Cell number 'ncl' points (i1,i2,i3) cell.
      icel(0,ncl)=0
      icel(1,ncl)=i1
      icel(2,ncl)=i2
      icel(3,ncl)=i3
      do 10 i=1,3
   10 a(i,ncl)=dble(i1)*r(i,1)+dble(i2)*r(i,2)+dble(i3)*r(i,3)
c     --- Count the total number of pairs for which distance
c         between the pair atoms is shorter than 'dlimit'
      dlimp=dlimit**2+small
      ndata=0
      do 20 ic=1,ncl
      do 20 ia1=1,natm
      do 20 ia2=ia1,natm
c     --- If the partner atom is located in a different cell and also
c         is not sitting on a crystallographically equivalent site,
c         the pair obtained by exchanging these atoms also should
c         be listed up. However, if the pair is within the
c         same cell, it is listed up only once.
      incr=2
      if(ic .eq. 1 .or. ia1 .eq. ia2) incr=1
      dst=0d0
      do 30 i=1,3
   30 dst=dst+(atmicp(i,ia2)-atmicp(i,ia1)+a(i,ic))**2
   20 if(dst .lt. dlimp .and. dst .gt. small) ndata=ndata+incr
c     --- Get the inverse mapping of 'irotat'.
      call clrari(invrot,natm*48)
      do 40 i=1,natm
   40 invrot(i,1)=i
      do 50 iop=2,48
      if(isymop(iop) .eq. 1) then
      do 60 i=1,natm
   60 invrot(irotat(i,iop),iop)=i
      endif
   50 continue
c
c     --- The main procedure starts here.
      call clrari(ipair,4*npairx)
      npair=0
      mjtbl=0
      ncel=0
      do 70 it1=1,ntyp
      call clrari(lookup,(2*nc(1)+1)*(2*nc(2)+1)*(2*nc(3)+1)*natm)
c     --- look for an atom of specified type.
      ia1=iatm(it1)
c     --- only one atom ia1 of type it1 is considered
      do 70 nl=1,ncl
      i1=icel(1,nl)
      i2=icel(2,nl)
      i3=icel(3,nl)
      do 70 it2=it1,ntyp
c     --- on the other hand, calculate all atoms ia2 of type=it2
      do 70 ia2=1,natm
      if(itype(ia2) .eq. it2) then
      dst=0d0
      do 80 i=1,3
      t(i,1)=a(i,nl)+atmicp(i,ia2)-atmicp(i,ia1)
      t(i,2)=t(i,1)
   80 dst=dst+t(i,1)**2
      if(lookup(i1,i2,i3,ia2) .eq. 0 .and. dst .gt. small .and.
     &       dst .lt. dlimp) then
c     --- a new partner found
c     --- list up all the components that are on the same site
c
      if(icel(0,nl) .eq. 0) then
c     --- This cell has not yet been registered and must be marked
c         as "to be calculated". Increment the number 'ncel' and
c         give an index 'ncel' to this cell.
      ncel=ncel+1
      icel(0,nl)=ncel
      if(ncel .gt. ncelx) call errtrp(1,'apairs','ncel too large')
      do 90 i=1,3
   90 d(i,ncel)=a(i,nl)
      endif
c     --- cartesian coordinates of the cell are given by 'd'.
c
c     write(*,'(a)')'partner atoms'
c     write(*,'(a,3i3,a,2i2)')'cell=',i1,i2,i3,' type,atom=',it2,ia2
      npair=npair+1
      if(npair .gt. npairx) call errtrp(1,'apairs','npair too large')
c     --- the first atom of the pair is 'ia1' sitting at the cell 1,
c         i.e. (0,0,0), and the second atom is 'ia2' at the cell 'nl'
      ipair(1,npair)=ia1
      ipair(2,npair)=ia2
      ipair(3,npair)=icel(0,nl)
      ipair(4,npair)=0
c     --- dval is used for sorting with respect to it1, it2, and dst
      dval(npair)=dble(ntyp*(it1-1)+it2)+sqrt(dst)/dlimit
c     --- list up all the equivalent pairs.
      ia3=0
      if(it1 .eq. it2) then
      ichk=2
c     --- As long as (1,2) pair is regarded as equivalent to (2,1) pair,
c         as is the case where the anti-symmetric exchange interaction
c         is neglected, it is reasonable to first find out the partner
c         that forms such a pair, namely find x for which (1,x) is
c         equivalent to (2,1). This procedure is necessary because
c         some of those pairs cannot be detected by symmetry
c         consideration alone. In most cases x=2, but in some cases
c         x=2 does not satisfy the above condition.
      if(ia1 .eq. ia2) then
c     --- The two atoms has the same position at the different cells.
c         The direction connecting two atoms must have been reversed.
      ia3=ia2
      do 100 i=1,3
  100 t(i,2)=-t(i,1)
      else
      do 110 iop=2,48
      if(irotat(ia1,iop) .eq. ia2) then
      call rotatm(t(1,1),t(1,2),protat(1,iop),1)
      do 120 i=1,natm
  120 if(irotat(i,iop) .eq. ia1) ia3=i
c     --- The direction connecting two atoms must have been reversed.
      do 130 i=1,3
  130 t(i,2)=-t(i,2)
      go to 140
      endif
  110 continue
      call errtrp(3,'apairs','symmetry not consistent with irotat')
  140 continue
c     --- Now it is found that the (ia1,ia2) pair is equivalent to
c         the (ia1,ia3) pair, which may not be detected using the
c         symmetry properties alone.
      endif
      else
      ichk=1
      endif
c     --- Now, mark up all the equivalents pairs as 'already done'.
      do 150 iop=1,48
      if(irotat(ia1,iop) .eq. ia1) then
c     --- This position has the symmetry 'iop'
      call rotatm(t,s,protat(1,iop),ichk)
c     --- It is neccessary to check both (ia1,ia2) and (ia1,ia3)
c         pairs because those two are not connected by symmetry.
c         In the case that ia3 does not exist, i.e. it1 and it2 are
c         different, ia3 has been set 0.
      do 160 i=1,natm
c     --- pairs equivalent to the pair (ia1, ia2)
      if(irotat(i,iop) .eq. ia2) itgt(1)=i
c     --- pairs equivalent to the pair (ia1, ia3)
  160 if(irotat(i,iop) .eq. ia3) itgt(2)=i
      do 170 j=1,ichk
c     --- To be on the safe side, the atomic position in the unit cell
c         is biased during the process finding out the cell where the
c         target atom resides. s(i,j), after being biased, becomes the
c         vector s2 connecting two cells.
      do 180 i=1,3
  180 s2(i)=atmicp(i,ia1)+s(i,j)-atmicp(i,itgt(j))
c     --- Find the cell where the rotated position resides, i.e. find
c         'icn' that satisfies s2=icn(1)*r(1)+icn(2)*r(2)+icn(3)*r(3),
c         where r are the primitive translation vectors.
      do 190 i=1,3
      cn(i)=0d0
      do 190 n=1,3
  190 cn(i)=cn(i)+s2(n)*g(n,i)
      do 200 i=1,3
  200 icn(i)=floor(cn(i)+1d-3)
c     --- The targeted atom must be 'itgt'.
c         Mark up this position as 'already done'.
      if(abs(icn(1)) .le. nc(1) .and. abs(icn(2)) .le. nc(2) .and.
     &   abs(icn(3)) .le. nc(3)) then
      if(lookup(icn(1),icn(2),icn(3),itgt(j)) .eq. 0) then
      lookup(icn(1),icn(2),icn(3),itgt(j))=1
      ipair(4,npair)=ipair(4,npair)+1
c     --- Partner atom is atom 'ia2' in the cell 'nl'. J_ij of this
c         pair is labeled as 'npair'
c
c     --- 'jtbl' lists up all the pairs that are crystallographically
c         equivalent to (ia1,ia2) pair.
c     --- At this stage the program has searched only the
c         pairs whose first atom is given by iatm(ityp).
c         In order to tabulate all the pairs whose first atom is
c         within the first cell (the central cell), it is necessary to
c         further find out the pairs whose first atom is not restricted
c         to iatm(ityp) but yet equivalent to the pair that has
c         iatm(ityp) as one of paired atoms.
      do 210 ja1=1,natm
c     --- The followintg procedure applied on the atom 'ja1' that
c         is in the first cell and of the same type of atom ia1.
c         This procedure also includes the case where ja1 =ia1,
c         i.e., the case that has been already found and designated
c         as the representative of pairs of the same kind.
      if(itype(ja1) .eq. it1) then
c     --- Find a symmetry operation that maps 'ia1' to 'ja1'
      do 220 iop2=1,48
      if(invrot(ia1,iop2) .eq. ja1) go to 230
  220 continue
      write(*,*)'ja1,itype(ja1),it1',ja1,itype(ja1),it1
      call errtrp(1,'apairs','target atom not found')
c     --- The second atom 'itgt(j)' of the pairs will be mapped to
c         ja2 by the same symmetry operation.
  230 ja2=invrot(itgt(j),iop2)
c     --- Find the cell where the mapped atom ja2 is sitting. This
c         cannot be done by rotating the cell where the atom ia2 is
c         residing becuase the rotation could associate a translation.
      call rotatm(s(1,j),s2,protat(1,iop2),1)
      do 240 i=1,3
  240 s2(i)=s2(i)+atmicp(i,ja1)-atmicp(i,ja2)
      do 250 i=1,3
      cn(i)=0d0
      do 250 n=1,3
  250 cn(i)=cn(i)+s2(n)*g(n,i)
      do 260 i=1,3
  260 icn(i)=floor(cn(i)+1d-3)
      if(abs(icn(1)) .le. nc(1) .and. abs(icn(2)) .le. nc(2) .and.
     &   abs(icn(3)) .le. nc(3)) then
c     --- Find the cell number 'nclt' that corresponds to the
c         (icn(1),icn(2),icn(3)) cell.
      do 270 i=1,3
  270 icl(i)=mod(icn(i)+2*nc(i)+1,2*nc(i)+1)
      nclt=(icl(1)*(2*nc(2)+1)+icl(2))*(2*nc(3)+1)+icl(3)+1
c     --- Also find the cell 'inclt' that is on the opposite direction.
      do 280 i=1,3
  280 icl(i)=mod(-icn(i)+2*nc(i)+1,2*nc(i)+1)
      inclt=(icl(1)*(2*nc(2)+1)+icl(2))*(2*nc(3)+1)+icl(3)+1
c     --- If the both atoms are in the same unit cell (cell-1)
c         and are the same type, (ja1,ja2) with ja1>ja2 is excluded
c         since (ja1,ja2) is equivalent to (ja2,ja1). If these
c         atoms belong to different types, only it2>it1 cases
c         are accepted (revrsed pairs are treated in the next).
      if(nclt .ne. 1 .or. (ja2 .gt. ja1 .and. it1 .eq. it2)
     &   .or. it2 .gt. it1) then
      mjtbl=mjtbl+1
      if(mjtbl .gt. mjtblx) call errtrp(1,'apairs','mjtbl too large')
      jtbl(1,mjtbl)=ja1
      jtbl(2,mjtbl)=ja2
      jtbl(3,mjtbl)=nclt
      jtbl(4,mjtbl)=npair
c     --- If the two atoms ja1 and ja2 are sitting on the different
c         cells, it is needed to list up the combination where
c         ja1 and ja2 are exchanged, putting ja2 on the first cell.
      if(inclt .ne. 1 .and. it2 .gt. it1) then
      mjtbl=mjtbl+1
      if(mjtbl .gt. mjtblx) call errtrp(1,'apairs','mjtbl too large')
      jtbl(1,mjtbl)=ja2
      jtbl(2,mjtbl)=ja1
      jtbl(3,mjtbl)=inclt
      jtbl(4,mjtbl)=npair
      endif
      endif
      endif
      endif
  210 continue
      endif
      endif
  170 continue
      endif
  150 continue
      endif
      endif
   70 continue
c     --- If the procedure be correctly done one might expect
c         ithat mjtbl=ndata must be satisfied. However, this is not
c         correct when the nc's are set such that some inter-cell
c         interactions are dropped. This is somehow dangerous because
c         the missing inter-cell interaction may be caused by atoms
c         located at very short distance. To avoid such a situation
c         if any nc's are set zero, the atomic geometry should be
c         carefullly checked. Also something wrong might be detected
c         by difference in mjtbl and ndata.
      write(*,'(a,i8,a,i6,a,i6)')
     &  '   npair=',npair,'  mjtbl=',mjtbl,'  ndata=',ndata
      if(mjtbl .ne. ndata)
     &  call errtrp(2,'apairs','values of mjtbl and ndata different')
c
c     --- Sort the data in the ascending order of dval(it1,it2,dst).
      do 290 i=1,npair
  290 map(i)=i
      do 300 j=2,npair
      swp=dval(j)
      mswp=map(j)
      do 310 k=1,4
  310 iswp(k)=ipair(k,j)
      do 320 i=j-1,1,-1
      if(dval(i) .lt. swp+1d-10) go to 330
      dval(i+1)=dval(i)
      map(i+1)=map(i)
      do 320 k=1,4
  320 ipair(k,i+1)=ipair(k,i)
      i=0
  330 dval(i+1)=swp
      map(i+1)=mswp
      do 300 k=1,4
  300 ipair(k,i+1)=iswp(k)
c     --- 'imap' gives the inverse mapping of 'map'
      do 340 i=1,npair
  340 imap(map(i))=i
c     --- Convert the index that refers the representative
c         pair for which J_ij is actually to be calculated.
      do 350 i=1,mjtbl
  350 jtbl(4,i)=imap(jtbl(4,i))
c
c     --- various optional output 1
c     do 430 i=1,npair
c     dst=(dval(i)-dble((itype(ipair(1,i))-1)*ncmpx
c    &         +itype(ipair(2,i))))*dlimit
c 430 write(*,'(1x,5i5,f10.5,2x,3f10.5)')
c    &  i,(ipair(j,i),j=1,4),dst,(d(j,ipair(3,i)),j=1,3)
      write(*,'(a,i4,a,i8)')'   ncel=',ncel,'  npair=',npair
cc    do 440 i=1,mjtbl
c     do 440 i=1,100
c 440 write(*,'(a,i3,a,i3,a,i3,a,i4)')'  atom1=',jtbl(1,i)
c    &  ,'  atom2=',jtbl(2,i),'  cell=',jtbl(3,i),'  indx=',jtbl(4,i)
c
c     --- optional output 2
c     write(*,'(//t6,a,t14,a,t20,a,t33,a,t44,a)')
c    & 'pair','atom1','atom2','cell','index'
c     do 400 i=1,mjtbl
c     ic=jtbl(3,i)
c 400 write(*,'(2x,i6,2x,2i6,6x,a,3i3,a,2x,i4)')
c    &  i,jtbl(1,i),jtbl(2,i),'(',(icel(j,jtbl(3,i)),j=1,3)
c    & ,' )',jtbl(4,i)
c     do 450 i=1,mjtbl
c     j=jtbl(4,i)
c     dst=(dval(j)-dble((itype(ipair(1,j))-1)*ncmpx
c    &         +itype(ipair(2,j))))*dlimit
c 450 write(*,'(a,3i4,f12.7)')'ia1,ia2,ic,dist',jtbl(1,i),jtbl(2,i),
c    &   jtbl(3,i),dst
      deallocate(lookup,dval,map,imap,a,aa,invrot,stat=ierr)
      end
