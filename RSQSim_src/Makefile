CMD	= runRSQSim
CC	= gcc
#CC      = mpicc
FC      = gfortran

#CSRCS   = $(wildcard *.c) 
CSRCS   = runRSQSim.c RSQSim.c allocFaultModel.c allocKRDOkadaFaultPatch.c \
          approxTgf.c writeSnapshots.c \
          calcStiffnessMatrix.c d1.c darray2.c dc3d.c dc3d0.c ddot1.c deformationToStress.c \
          farray2.c findCorners.c findNormalAndSlip.c fixedPoint.c getDunavantPtsWeights.c \
          initFaultModel.c \
          nlines.c packFaultModel.c packParams.c projectStress.c readFault.c \
          setDdot3.c strflds.c theta_0.c theta_1.c timeStep.c transitionTimes.c \
          transitions.c xfopen.c xmalloc.c zbrent.c utils_x.c fAfunc.c
ifeq (,$(findstring getParams.c,$(CSRCS)))
CSRCS += getParams.c
endif
#FSRCS = $(wildcard *.f)
FSRCS   = tgf_merged_s07.f
SRCS    = $(CSRCS) $(FSRCS)
PPRSRCS   = $(wildcard RSQSimPostProcess/*/*)
BFMRSRCS   = $(wildcard RSQSimBuildFaultModels/*/*)
#HDRS    = $(wildcard *.h)
HDRS    = KRDOkada.h RSQSim.h RSQSim_dataStructures.h RSQSim_defaultParams.h \
          RSQSim_functionDec.h RSQSim_params.h dc3d.h xfile.h xmalloc.h
COBJS	= $(patsubst %.c, %.o, $(CSRCS))
FOBJS	= $(patsubst %.f, %.o, $(FSRCS))
OBJS    = $(COBJS) $(FOBJS) 
DBXFLAG = -O
CFLAGS	= -Wall $(DBXFLAG)  -m64  #  -DUSE_SERIAL_IO
FFLAGS = $(DBXFLAG) -m64
CLIBS   = -L/sw/lib -lm -lgfortran
#RARCH   = --arch x86_64

all: 
	make RSQSimPostProcessT
	make RSQSimBuildFaultModelsT
	make $(CMD)
	make .RprofileT
	make singlePlanarFault
	make singleStrikeSlip
	-make -C Utils all

getParams.c: $(OBJS) $(HDRS) Makefile getParams.c.in
	sed "s/@svnversion.*@/@svnversion: `svnversion`@/g" getParams.c.in > getParams.c

.RprofileT:
	sed "s@REPLACE@`pwd`@" .Rprofile.template > .Rprofile

$(OBJS): $(HDRS) 

$(CMD): $(OBJS) Makefile $(HDRS)
	$(CC) $(CFLAGS) $(OBJS) $(CLIBS) -o $@

RSQSimPostProcessT: $(PPRSRCS)
	-rm .Rprofile
#	R --vanilla -e 'setwd("RSQSimPostProcess/R");library(devtools);document()'
	R $(RARCH) --vanilla CMD build RSQSimPostProcess
	R $(RARCH) --vanilla CMD INSTALL --library=./Rlibs RSQSimPostProcess_1.0.tar.gz
	make .RprofileT

RSQSimBuildFaultModelsT: $(BFMRSRCS)
	-rm .Rprofile
	R --vanilla -e 'setwd("RSQSimBuildFaultModels/R");library(devtools);document()'
	R $(RARCH) --vanilla CMD INSTALL --library=./Rlibs RSQSimBuildFaultModels
	make .RprofileT

test:
	-rm TestCase/singleSS/eqs.test1.out
	cd TestCase/singleSS; ../../runRSQSim test1.in
	grep -v outFnameInfix TestCase/singleSS/eqs.test.out | grep -v ProgramVersion | grep -v StartTime | grep -v Hostname | grep -v CommandLine | grep -v RunDirectory > TestCase/singleSS/eqs.test.out.tmp
	grep -v outFnameInfix TestCase/singleSS/eqs.test1.out | grep -v ProgramVersion | grep -v StartTime | grep -v Hostname | grep -v CommandLine | grep -v RunDirectory > TestCase/singleSS/eqs.test1.out.tmp
	diff TestCase/singleSS/eqs.test.out.tmp TestCase/singleSS/eqs.test1.out.tmp
	rm TestCase/singleSS/eqs.test.out.tmp TestCase/singleSS/eqs.test1.out.tmp 

redoTestCase:
	-rm TestCase/singleSS/eqs.test.out
	cd TestCase/singleSS; ../../runRSQSim test.in
        
singlePlanarFault: BuildFaultModels/singlePlanarFault.c
	gcc -o BuildFaultModels/singlePlanarFault BuildFaultModels/singlePlanarFault.c -lm

singleStrikeSlip: BuildFaultModels/singleStrikeSlip.c
	gcc -o BuildFaultModels/singleStrikeSlip BuildFaultModels/singleStrikeSlip.c -lm

#tarball: 
#	-rm hypoDD.C.tgz
#	tar cvf - --dereference *.c *.h Makefile test | gzip > hypoDD.C.tgz
#
#tarball_notest: 
#	-rm hypoDD.C.tgz
#	tar cvf - --dereference *.c *.h Makefile | gzip > hypoDD.C.tgz

clean:
	-rm -f *.o *.mod

distclean:
	-rm -f $(CMD) $(OBJS) *.o *.mod core a.out *.lint *.out *.aux .Rprofile *.dvi *.log *.bbl *.blg *.ps getParams.c
	-rm -rf Rlibs/* RSQSimPostProcess/src-* RSQSimPostProcess_1.0.tar.gz
	-rm -f TestCase/*test1*out* TestCase/*List
	-rm BuildFaultModels/singlePlanarFault BuildFaultModels/singleStrikeSlip
	make -C Utils clean
	-rm TestCase/TIMES_PARTITION.out TestCase/Screen.Time.out
	-rm TestCase/times_saved.out
	-rm RSQSimBuildFaultModels/man/*
