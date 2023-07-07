      subroutine wscell(r,atmicp,iatm,itype,natm,ntyp,ws,g)
c-----------------------------------------------------------------------
c     Given atomic positions, this program constructs the Wigner-Seitz
c     (Voronoi) cells for each type of atom and returns the radii
c     corresponding to the volumes of cells.
c     Currently the variable itype is not used. However, this may be
c     used to provide the specific feature of each atom.
c     Coded by H. Akai, 29 May. 2019, Tokyo.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(nz=14)
c     --- npx=60 is probably enough. nv could be as large as 40:
c         nvx=40 or larger should be specified.
      real*8 atmicp(3,natm),ws(ntyp),r(3,3),a(3,3),b(3),tmp(4),g(3,3)
      real*8,allocatable::p(:,:),d(:),vx(:,:),y(:,:)
      integer iatm(ntyp),itype(natm),itmp(2),indx(3)
      integer,allocatable::iface(:,:),lmn(:,:)
      data tiny/1d-8/,dmrgn/2.2d0/
      pi=4d0*atan(1d0)
      vc=(r(2,1)*r(3,2)-r(3,1)*r(2,2))*r(1,3)
     &  +(r(3,1)*r(1,2)-r(1,1)*r(3,2))*r(2,3)
     &  +(r(1,1)*r(2,2)-r(2,1)*r(1,2))*r(3,3)
      dlimit=dmrgn*(7.5d-1*vc/pi/dble(natm))**(1d0/3d0)
      do 10 it=1,ntyp
      write(*,*)
      write(*,*)'--- type=',it,' ---'
      ws(it)=0d0
      ia=iatm(it)
      np=0
      do 20 l=-2,2
      do 20 m=-2,2
      do 20 n=-2,2
      do 20 ib=1,natm
      do 30 i=1,3
   30 b(i)=dble(l)*r(i,1)+dble(m)*r(i,2)+dble(n)*r(i,3)
     &       +atmicp(i,ib)-atmicp(i,ia)
      dd=5d-1*sqrt(b(1)**2+b(2)**2+b(3)**2)
   20 if(dd .gt. tiny .and. dd .lt. dlimit) np=np+1
c
c     --- the values of nvx and nfx cannot be fixed in advance.
c         the following values might or might not be safe.
      nvx=2*np
      nfx=16*np
c
      allocate(p(3,np),d(np),vx(3,nvx),iface(2,nfx),y(4,nz)
     &        ,stat=ierr)
      if(ierr .ne. 0) call errtrp(1,'wscell','allocation fails')
      np=0
      do 40 l=-2,2
      do 40 m=-2,2
      do 40 n=-2,2
      do 40 ib=1,natm
      do 50 i=1,3
   50 b(i)=dble(l)*r(i,1)+dble(m)*r(i,2)+dble(n)*r(i,3)
     &       +atmicp(i,ib)-atmicp(i,ia)
      dd=5d-1*sqrt(b(1)**2+b(2)**2+b(3)**2)
      if(dd .gt. tiny .and. dd .lt. dlimit) then
      np=np+1
      do 60 i=1,3
   60 p(i,np)=b(i)
      d(np)=2d0*dd**2
      endif
   40 continue
      nf=0
      nv=0
      do 70 l=1,np-2
      do 70 m=l+1,np-1
      do 70 n=m+1,np
c     --- three planes are chosen.
      do 80 i=1,3
      a(1,i)=p(i,l)
      a(2,i)=p(i,m)
   80 a(3,i)=p(i,n)
c     --- get position of the vertex formed by the three planes
      b(1)=d(l)
      b(2)=d(m)
      b(3)=d(n)
      call rludcm(a,3,indx,dd)
      if(abs(dd*a(1,1)*a(2,2)*a(3,3)) .lt. tiny) go to 70
      call rlubks(a,3,indx,b)
c     --- check if the vertex and (0,0,0) is in the same side of all
c         other planes.
      do 90 j=1,np
      if((j-l)*(j-m)*(j-n) .ne. 0 .and.
     &   p(1,j)*b(1)+p(2,j)*b(2)+p(3,j)*b(3)-d(j) .gt. tiny) go to 70
   90 continue
c     --- the vertex and (0,0,0) is in the same side of all
c         other planes.
      do 100 j=1,nv
c     --- sometimes more than three planes meet at a single vertex.
      iv=j
      if((b(1)-vx(1,j))**2+(b(2)-vx(2,j))**2+(b(3)-vx(3,j))**2
     &   .lt. tiny) go to 110
  100 continue
c     --- this vertex must be a new one.
      nv=nv+1
      iv=nv
      if(nv .gt. nvx) call errtrp(1,'wscell','nv too large')
c     --- now a vertex that is given by three planes (faces) (l,m,n)
c         is fixed.
      do 120 i=1,3
  120 vx(i,nv)=b(i)
c     --- store three faces, l, m, and n, and the corresponding vertex.
  110 if(nf+3 .gt. nfx) call errtrp(1,'wscell','nf too large')
      iface(1,nf+1)=l
      iface(1,nf+2)=m
      iface(1,nf+3)=n
c     --- iface(2,*) contains the vertex associated to the faces.
      do 130 j=1,3
  130 iface(2,nf+j)=iv
      nf=nf+3
   70 continue
c     --- now all the vertices are found.
c     write(*,'(a,i3)')'number of planes considerd =',np
c     write(*,'(a,i3)')'number of faces forming vertices=',nf
c     write(*,'(a,i3)')'number of vrtices =',nv
c     write(*,'(i3,3e13.5)')(n,(vx(i,n),i=1,3),n=1,nv)
c
c     --- first, sort iface(*,face number) in the ascending order of the
c         vertex number.
      do 140 j=2,nf
      do 150 l=1,2
  150 itmp(l)=iface(l,j)
      do 160 i=j-1,1,-1
      if(iface(2,i) .le. itmp(2)) go to 170
      do 180 l=1,2
  180 iface(l,i+1)=iface(l,i)
  160 continue
      i=0
  170 do 190 l=1,2
  190 iface(l,i+1)=itmp(l)
  140 continue
c     write(*,'(3i3)')(i,(iface(l,i),l=1,2),i=1,nf)
c
c     --- then, sort iface(*,face number) in the ascending order of the
c         face number.
      do 200 j=2,nf
      do 210 l=1,2
  210 itmp(l)=iface(l,j)
      do 220 i=j-1,1,-1
      if(iface(1,i) .le. itmp(1)) go to 230
      do 240 l=1,2
  240 iface(l,i+1)=iface(l,i)
  220 continue
      i=0
  230 do 250 l=1,2
  250 iface(l,i+1)=itmp(l)
  200 continue
c     write(*,'(3i3)')(i,(iface(l,i),l=1,2),i=1,nf)
c     --- discard the equivalent (vertex, face) pairs.
      nn=1
      do 260 i=1,nf-1
      if(iface(1,i) .ne. iface(1,i+1) .or.
     &   iface(2,i) .ne. iface(2,i+1)) then
      nn=nn+1
      do 270 l=1,2
  270 iface(l,nn)=iface(l,i+1)
      endif
  260 continue
c     --- now number of data in iface is reduced from nf to nn
      nf=nn
c     write(*,'(3i3)')(i,(iface(l,i),l=1,2),i=1,nf)
c     --- istrt points the start of a new face.
      istrt=1
c     --- repeat the following procedure until all the faces are
c         processed.
      do 280 irpt=1,nf-1
c     --- check if the subsequent iface points the same face
c         as the current one.
      nx=1
      do 290 i=1,nf-istrt
c     write(*,*)i,nx,iface(1,istrt+i),iface(1,istrt)
      if(iface(1,istrt+i) .ne. iface(1,istrt)) exit
c     --- nx should be determined here to ensure the
c         correct number of vertex irrespective of
c         the above inequality fulfilled or not.
  290 nx=i+1
c     write(*,'(a,i3,a,i3)')'istrt=',istrt,' nx for this face=',nx
      if(nx-1 .gt. nz) call errtrp(1,'wscell','nz too small')
c     --- if mx<3, no face can be formed. this corresponds
c         to the case where the plane has a single vertex or an edge.
      if(nx .ge. 3) then
c     --- nx vertices are associated to this face. Some
c         verticies may coincide with each other's.
c     --- y's are nx-1 vectors connecting the base vertex and other
c         vertices,
      ibase=iface(2,istrt)
      do 300 l=1,3
  300 a(l,1)=vx(l,ibase)
      itrgt=iface(2,istrt+1)
      do 310 l=1,3
  310 y(l,1)=vx(l,itrgt)-a(l,1)
      yb=sqrt(y(1,1)**2+y(2,1)**2+y(3,1)**2)
      y(4,1)=0d0
      do 320 i=2,nx-1
      itrgt=iface(2,istrt+i)
      do 330 l=1,3
  330 y(l,i)=vx(l,itrgt)-a(l,1)
      ya=sqrt(y(1,i)**2+y(2,i)**2+y(3,i)**2)
c     --- construct a vector that gives the direction of a plane
c         spanned by two vectors y_1 and y_i.
      b(1)=y(2,1)*y(3,i)-y(3,1)*y(2,i)
      b(2)=y(3,1)*y(1,i)-y(1,1)*y(3,i)
      b(3)=y(1,1)*y(2,i)-y(2,1)*y(1,i)
c     --- obtain the sign of direction of the directional vector.
      sgn=sign(1d0,b(1)*a(1,1)+b(2)*a(2,1)+b(3)*a(3,1))
c     --- y(4,i) contains the angle between two vectors y_i and y_1.
c     --- the angle has a sign, depending on whether y_i locates to the
c         left or to the right of y_1.
      y(4,i)=acos(min(1d0,(ya**2+yb**2
     &     -(y(1,i)-y(1,1))**2-(y(2,i)-y(2,1))**2-(y(3,i)-y(3,1))**2)
     &     /(2d0*ya*yb)))*sgn
c     write(*,'(a,i3,1p5e13.6)')'i,ya,yb,sgn,y',i,ya,yb,sgn,y(4,i)
c     write(*,'(a,2i3,4f120.5)')'i,itrgt,y(4,i)',i,itrgt,y(4,i)
  320 continue
c    --- sort y(*,i) in ascending order of the angle y(4,*).
      do 340 j=2,nx-1
      do 350 l=1,4
  350 tmp(l)=y(l,j)
      do 360 i=j-1,1,-1
      if(y(4,i) .lt. tmp(4)) go to 370
      do 380 l=1,4
  380 y(l,i+1)=y(l,i)
  360 continue
      i=0
  370 do 390 l=1,4
  390 y(l,i+1)=tmp(l)
  340 continue
c     write(*,*)'number of target vertices=',nx-1
c     write(*,'(i3,4f12.5)')(i,(y(l,i),l=1,4),i=1,nx-1)
c
c     --- calculate the volume of parallelpiped spanned by three
c         vectors vx(ibase), vx(ibase)+y(i), and vx(ibase)+y(i+1)
c         for i=1,nx-2
      do 400 i=1,nx-2
      do 410 l=1,3
      a(l,2)=a(l,1)+y(l,i)
  410 a(l,3)=a(l,1)+y(l,i+1)
  400 ws(it)=ws(it)
     &     +abs(a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
     &         +a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3))
     &         +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3)))/6d0
      endif
c     --- The new face will start from new istrt=istrt+nx
      istrt=istrt+nx
      if(istrt .gt. nf-1) exit
  280 continue
      write(*,'(a,i3)')'number of vertices =',nv
      do 440 n=1,nv
      dvrt=sqrt(vx(1,n)**2+vx(2,n)**2+vx(3,n)**2)
      do 450 i=1,3
  450 vx(i,n)=vx(i,n)+atmicp(i,ia)
      do 460 i=1,3
      b(i)=0d0
      do 470 j=1,3
  470 b(i)=b(i)+vx(j,n)*g(j,i)
c 460 b(i)=b(i)-floor(b(i))
  460 continue
c 440 if(abs(dvrt-0.31115) .lt. 1d-3)
c 440 if(abs(dvrt-0.43301) .lt. 1d-3)
c    & write(*,'(i3,f12.5,3(f12.7,a1))')
  440 write(*,'(i3,f12.5,3(f12.7,a1))')
     &    n,dvrt,b(1),'a',b(2),'b',b(3),'c'
      deallocate(p,d,vx,iface,y,stat=ierr)
   10 continue
      vol=0d0
      do 420 i=1,natm
  420 vol=vol+ws(itype(i))
      write(*,*)vol,vc
      if(abs(vol-vc) .gt. tiny)
     &  call errtrp(1,'wscell','procedure fails')
       do 430 i=1,ntyp
  430 ws(i)=(7.5d-1*ws(i)/pi)**(1d0/3d0)
      return
      end
