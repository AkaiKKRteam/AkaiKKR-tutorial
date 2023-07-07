#include <unistd.h>
#include <string.h>
#include <sys/utsname.h>
#include <stdio.h>
int numcor_(char *os,char *host,char *arch)
/*
------------------------------------------------------------------------
   This program return the number of  processors
   that can be used.
   For some systems the following values are set:
   _SC_NPROCESSORS_CONF = 57
   _SC_NPROCESSORS_ONLN = 58
   However, these values could vary from system to system.
   dsp controls the optional output as for the operating system.
   Coded by H. Akai, Tokyo, Jul 5, 2022
------------------------------------------------------------------------
*/
{
   struct utsname sys;
   int numc,dsp=0;
   numc=sysconf(_SC_NPROCESSORS_ONLN);
   if(uname(&sys) == 0) {
     strcpy(os,sys.sysname);
     strcpy(host,sys.nodename);
     strcpy(arch,sys.machine);
     if(dsp) {
       printf(" OS: %s\n",os);
       printf(" Host: %s\n",host);
       printf(" Hardware: %s\n\n",arch);}
/* For fx700 or fx1000 2 nodes should be reserved for the system */
     if(strstr(sys.machine,"aarch64") != NULL)
       {numc=numc-2;}
    }
   return numc;
/*
   for availabe number of processors
   printf("%d\n",_SC_NPROCESSORS_ONLN);
   Use _SC_NPROCESSORS_CONF for maximum number of preseccors
*/
}
