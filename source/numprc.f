      function numprc() 
c-----------------------------------------------------------------------
c     This function return the available number of cores.
c     Codec by H. Akai, Tokyo, 5 July 2022
c-----------------------------------------------------------------------
      interface
      function get_nprocs() bind(c)
      integer:: get_nprocs
      end function
      end interface
      numprc=get_nprocs()
      end
