CC      = /usr/local/gfortran/bin/gcc
FC      = /usr/local/gfortran/bin/gfortran

DBXFLAG = -O
CFLAGS	= -Wall $(DBXFLAG) -DUSE_MPI
FFLAGS = -Wall $(DBXFLAG) -arch x86_64
CLIBS   = -L/sw/lib -lm -lgfortran -lquadmath

all: 
ifneq ($(wildcard $(CC)),)
#	make useOkadaQ
endif

$(OBJS): $(HDRS) 

$(CMD): $(OBJS) Makefile $(HDRS)
	$(CC) $(CFLAGS) $(OBJS) $(CLIBS) -o $@

# need quadmath for the following
#useOkadaQ: useOkadaQ.o ../dc3dq.o ../dc3dq.h
#	$(CC) $(CFLAGS) useOkadaQ.o ../dc3dq.o $(CLIBS) -o $@

CC = mpicc
FC = gfortran
CLIBS   = -L/sw/lib -lm -lgfortran 
CFLAGS	= -Wall $(DBXFLAG) -DUSE_MPI -m64
FFLAGS = -Wall $(DBXFLAG) -m64
COBJS = ../allocFaultModel.o ../readFault.o ../dc3d.o ../xmalloc.o ../xfopen.o \
        ../allocKRDOkadaFaultPatch.o ../darray2.o ../findCorners.o ../findNormalAndSlip.o \
        ../nlines.o ../approxTgf.o ../dc3d0.o ../deformationToStress.o ../getDunavantPtsWeights.o \
        ../tgf_merged_s07.o

sdotF: sdotF.o $(COBJS)
	$(CC) $(CFLAGS) sdotF.o $(COBJS) $(CLIBS) -o $@

clean:
	-rm -f useOkadaQ useOkadaQ.o sdotF sdotF.o core a.out *.lint *.out *.aux *.dvi *.log *.bbl *.blg *.ps
