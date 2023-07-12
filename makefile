# ThecCommand invoking the fortran compiler of the
# system should be specified by "fort= ".
# The corresponding compiler option also should be given 
# by "flag= ".
#
#.PHONY: all
#all: $(program) gpd spc
# If openMP is not supported two lines for omp should be
# replaced by 
#omp =
#nomp =
#
#######################################
# specify the fortran compiler
#fort = ifort
fort = gfortran
#######################################
ifeq ($(fort),ifort)
 	flag = -O2 -mcmodel=medium
#	flag = -O2 -mcmodel=medium -traceback
#	flag = -O2 -xcore-avx512 -mcmodel=medium
#	flag = -O2 -mcmodel=medium -check uninit -ftrapuv -traceback
#	flag = -mcmodel=medium -ftrapuv -traceback
#	flag = -O2 -mcmodel=medium -check all -traceback
	omp = -qopenmp
	nomp = -qopenmp-stubs
else
	flag = -O2
	omp = -fopenmp
	nomp =
endif
#
#######################################
# In the case of ifort 16.0.0 or later
#######################################
#fort = ifort
#flag = -O2 -mcmodel=medium
#flag = -O2 -xcore-avx512 -mcmodel=medium -ftrapuv -traceback
#flag = -O2 -xcore-avx512 -mcmodel=medium
#flag = -O2 -mcmodel=medium -mkl -ftrapuv -traceback
#flag = -O2 -mcmodel=medium -ftrapuv -traceback
#flag = -O2 -mcmodel=medium
#omp = -qopenmp
#omp = -qopenmp
#nomp = -qopenmp-stubs
#
#######################################
# In the case of ifort 15.x.x or ealier
#######################################
#fort = ifort
#flag = -O2 -mcmodel=medium
##flag = -mp1 -i-dynamic -mcmodel=medium
# use openmp for old version of ifort
#omp = -openmp
#nomp = -openmp-stubs
#
program = specx
objs = \
source/$(program).o \
source/spmain.o \
source/alphac.o \
source/apairs.o \
source/asymbl.o \
source/atmcor.o \
source/atmicv.o \
source/atmmas.o \
source/atmnum.o \
source/atmrot.o \
source/aveari.o \
source/avearr.o \
source/banden.o \
source/biomix.o \
source/bravai.o \
source/brintg.o \
source/brvsym.o \
source/brvsy2.o \
source/broydn.o \
source/bzmesh.o \
source/bzmrot.o \
source/bzmrtf.o \
source/bzmsmb.o \
source/bzmsmf.o \
source/bzvrtx.o \
source/cartes.o \
source/cemesh.o \
source/cemesr.o \
source/cg.o \
source/cgc.o \
source/cgntcd.o \
source/cgntcs.o \
source/cgnwt.o \
source/cgtabl.o \
source/chebpc.o \
source/chkcnf.o \
source/chklat.o \
source/chkprm.o \
source/chktyp.o \
source/chleng.o \
source/chrasa.o \
source/chrdnc.o \
source/cinvrn.o \
source/cinvrx.o \
source/ckmesh.o \
source/clrarc.o \
source/clrari.o \
source/clrarr.o \
source/cmdbsl.o \
source/cnsole.o \
source/corada.o \
source/corcnf.o \
source/corgga.o \
source/corlsd.o \
source/cpaitr.o \
source/cpolin.o \
source/cprinv.o \
source/cpshft.o \
source/crtint.o \
source/cryrho.o \
source/cstatc.o \
source/cstate.o \
source/cubmtr.o \
source/curie.o \
source/dawson.o \
source/diagrs.o \
source/diffn.o \
source/dplini.o \
source/drivtv.o \
source/drvlud.o \
source/drvmsh.o \
source/dsenum.o \
source/dspcpa.o \
source/dspgm.o \
source/equarc.o \
source/equari.o \
source/equarr.o \
source/equarx.o \
source/eqvlat.o \
source/eranlb.o \
source/erfc.o \
source/erranc.o \
source/errtrp.o \
source/etaopt.o \
source/excev.o \
source/excg91.o \
source/excgga.o \
source/exclmm.o \
source/excmjw.o \
source/excpbe.o \
source/excpym.o \
source/excpyv.o \
source/excvbh.o \
source/excvwn.o \
source/exgga.o \
source/extorg.o \
source/fczero.o \
source/fd3.o \
source/finitn.o \
source/fintgr.o \
source/fldf.o \
source/fxspin.o \
source/gcor.o \
source/gengpt.o \
source/genrpt.o \
source/gensm.o \
source/getcrt.o \
source/getdtb.o \
source/getfil.o \
source/getnum.o \
source/getorg.o \
source/getrtm.o \
source/gettmp.o \
source/getwsr.o \
source/gintgr.o \
source/gmtrns.o \
source/gnpset.o \
source/gntcds.o \
source/gntcs.o \
source/gpsort.o \
source/green.o \
source/gsdatp.o \
source/gtchst.o \
source/guesse.o \
source/guessz.o \
source/hexmtr.o \
source/hhforb.o \
source/hypera.o \
source/ibrava.o \
source/ifkey.o \
source/inijip.o \
source/inqbrv.o \
source/jkkr.o \
source/julday.o \
source/kkrsed.o \
source/laguer.o \
source/lftfil.o \
source/lindx.o \
source/lsqftr.o \
source/lstcgt.o \
source/lubksb.o \
source/ludcmp.o \
source/lwcase.o \
source/madlng.o \
source/mdbsl.o \
source/mdexp2.o \
source/mkemap.o \
source/mseque.o \
source/mshatm.o \
source/numcor.o \
source/nesbet.o \
source/neutrl.o \
source/nfqlty.o \
source/nzelem.o \
source/phasea.o \
source/phaseb.o \
source/phasef.o \
source/pltmrk.o \
source/poisna.o \
source/poisnb.o \
source/polin0.o \
source/polint.o \
source/potenv.o \
source/prdmtc.o \
source/prdmtr.o \
source/prdrtm.o \
source/prmvec.o \
source/prntcs.o \
source/prntpn.o \
source/pyexch.o \
source/qromo2.o \
source/qvolum.o \
source/radial.o \
source/radnrl.o \
source/ratin0.o \
source/redata.o \
source/readin.o \
source/readk.o \
source/realh.o \
source/reconf.o \
source/reduce.o \
source/relred.o \
source/rgconv.o \
source/rlubks.o \
source/rludcm.o \
source/rmesha.o \
source/rmserb.o \
source/rmserr.o \
source/rndkpt.o \
source/rnmrdt.o \
source/rnmrdx.o \
source/rotall.o \
source/rotatm.o \
source/rotave.o \
source/rrcnv.o \
source/rsymeq.o \
source/sbrnch.o \
source/sbtime.o \
source/setarc.o \
source/setari.o \
source/setarr.o \
source/setomp.o \
source/sintgr.o \
source/spckkr.o \
source/stfact.o \
source/srtrns.o \
source/subscr.o \
source/swpari.o \
source/swparr.o \
source/tchmta.o \
source/tchstr.o \
source/tfp.o \
source/tltcry.o \
source/totalw.o \
source/transp.o \
source/trmkey.o \
source/ty2ity.o \
source/uclock.o \
source/udate.o \
source/utimer.o \
source/uutrns.o \
source/uxcor.o \
source/viomix.o \
source/vrotat.o \
source/vwncor.o \
source/wrtspc.o \
source/wscell.o \
source/xtoken.o \
source/zroots.o

$(program): $(objs) 
	@$(fort) -o $@ $(flag) $(omp) $(objs) $(libs)

.f.o:
	$(fort) $(flag) $(omp) -o $(<:.f=.o) -c $<

gpd: source/gpd.f
	$(fort) $(flag) -o $@ $<

spc: source/spc.f
	$(fort) $(flag) -o $@ $<
