      function ifkey(key,keys)
c-----------------------------------------------------------------------
c     Given key and keys this program returns if string key
c     is involved in string keys.
c     When 'remove' is set true, an additional
c     operation will take place: If 'ifkey' is true, the corresponding key
c     is removed from keys.
c
c     *** in the presnet setting the key is not removed.
c     Coded by H. Akai, 23 Nov. 2007, revised Mar. 29, 2018.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character keys*(*),key*(*),dmys*180,dmy*180
      logical ifkey, remove
      data remove/.false./
      ifkey=.false.
      dmys=keys
      dmy=key
      m=lftfil(dmys)
      n=lftfil(dmy)
      i=index(dmys(1:m),dmy(1:n))
      if(i .gt. 0) then
      ifkey=.true. 
      if(remove) then
      keys(i:i+n-1)=' '
      m=lftfil(keys)
      endif
      endif
      end
