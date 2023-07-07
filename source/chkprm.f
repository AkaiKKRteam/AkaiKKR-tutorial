      subroutine chkprm(r)
c-----------------------------------------------------------------------
c     Given three vectors r(:,1),r(:,2),r(:,3), this program checks if
c     these vectors are primitive vectors of minimum length. If it
c     is not the case the program replaces the original vectors with
c     new vectors that have the minimum length.
c     Coded by H. Akai, Tokyo, May 4, 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter(mx=5,zero=1d-8)
      real*8 r(3,3),g(3,3),v(3,3)
      real*8,allocatable::d(:),p(:,:),cn(:,:)
      integer idone(3)
      integer,allocatable::indx(:)
      allocate(d(8*mx**3),p(3,8*mx**3),cn(3,8*mx**3),indx(8*mx**3))
c     --- Generate the lattice points for which the distance from the
c         origin is less than rmx.
      rmx=0d0
      do 10 i=1,3
   10 rmx=max(rmx,r(1,i)**2+r(2,i)**2+r(3,i)**2)
      nr=0
      do 20 l=-mx,mx
      do 30 i=1,3
   30 v(i,1)=dble(l)*r(i,1)
      do 20 m=-mx,mx
      do 40 i=1,3
   40 v(i,2)=v(i,1)+dble(m)*r(i,2)
      do 20 n=-mx,mx
      do 50 i=1,3
   50 v(i,3)=v(i,2)+dble(n)*r(i,3)
      dd=v(1,3)**2+v(2,3)**2+v(3,3)**2
      if(dd .lt. rmx+zero) then
      nr=nr+1
      d(nr)=dd
      do 60 i=1,3
   60 p(i,nr)=v(i,3)
      endif
   20 continue
c     --- Sort the vectors p in ascending order of their length.
      indx(1)=1
      do 70 i=2,nr
      indx(i)=i
      do 70 j=i-1,1,-1
      if(d(i) .lt. d(indx(j))) then
      indx(j+1)=indx(j)
      indx(j)=i
      endif
   70 continue
c     write(*,'(2i3,f12.7)')(i,indx(i),d(indx(i)),i=1,nr)
c     --- generate reciprocal lattice vectors g.
      call rgconv(r,g,vc)
      id=0
      do 80 n=2,nr
c     --- expand p(:,n) into r(:,i)
      do 90 i=1,3
      cn(i,n)=0d0
      do 90 j=1,3
   90 cn(i,n)=cn(i,n)+p(j,indx(n))*g(j,i)
c     write(*,'(i3,3f12.7)')n,(cn(i,n),i=1,3)
      do 100 i=1,3
      if(abs(cn(i,n)-1d0) .lt. zero .and.
     &   abs(cn(mod(i,3)+1,n)) .lt. zero .and.
     &   abs(cn(mod(i+1,3)+1,n)) .lt. zero) then
c     --- p(:,n) coincides one of r(:,i). New r(:,i) should
c         be the same as old r(:,i).
      id=id+1
      idone(id)=i
      exit
      endif
  100 continue
      if(id .eq. 3) exit
   80 continue
      if(id .lt. 3) then
      do 110 n=2,nr
      if(id .eq. 1) then
c     --- Is linearly independent of r(:,idone(1))?
c         Check if p x r_1 ist 0 or not.
      if((p(2,indx(n))*r(3,idone(1))-p(3,indx(n))*r(2,idone(1)))**2
     & +(p(3,indx(n))*r(1,idone(1))-p(1,indx(n))*r(3,idone(1)))**2
     & +(p(1,indx(n))*r(2,idone(1))-p(2,indx(n))*r(1,idone(1)))**2
     &   .lt. zero) go to 110
      else if(id .eq. 2) then
c     --- Is linearly independent of r(:,idone(1)) and r(:,idone(2))?
c         Check if (p x r_1) r_2 is 0 or not.
      if(abs((p(2,indx(n))*r(3,idone(1))-p(3,indx(n))*r(2,idone(1)))
     &  *r(1,idone(2))
     & +(p(3,indx(n))*r(1,idone(1))-p(1,indx(n))*r(3,idone(1)))
     &  *r(2,idone(2))
     & +(p(1,indx(n))*r(2,idone(1))-p(2,indx(n))*r(1,idone(1)))
     &  *r(3,idone(2))) .lt. zero) go to 110
      endif
c     --- New r found. Assign it as r(:,i) where i is one of
c         yet unassigned number among 1 through 3.
      if(id .eq. 0) ii=1
      if(id .eq. 1) ii=mod(idone(1),3)+1
      if(id .eq. 2) ii=6-idone(1)-idone(2)
      id=id+1
      idone(id)=ii
      do 120 i=1,3
  120 r(i,ii)=p(i,indx(n))
      if(id .eq. 3) exit
  110 continue
      endif
c     write(*,'(i3,3f12.7)')(i,(r(j,i),j=1,3),i=1,3)
      deallocate(d,p,cn,indx)
      end
