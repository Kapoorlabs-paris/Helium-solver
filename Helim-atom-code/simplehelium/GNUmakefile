PROG	 	= TDSE
PACKAGE  	= $(LIBDIR)/libqprop.a

ESQLBASICS    	= 
C++BASICS     	= TDSE.cc wavefunction.cc hamop.cc grid.cc fluid.cc
CBASICS      	= 
PROGBASICS      = TDSE

LOCALLIBS	=$(LIBDIR)/libqprop.a 

include         ../GNUmakefile.tmpl

#	User definable options

MYLIBS          = -lqprop -lm

COPTI           = -O8 #-pg
C++OPTI         = -O8 #-pg
C++SPECIAL	= -funroll-all-loops -finline-functions -fexpensive-optimizations -ffast-math -Wno-deprecated#-DTEST
LDFLAGS         = -Xlinker -defsym -Xlinker MAIN__=main\
                   # -L/opt/mpich/ch-p4/lib64/ #-pg #-gstabs+

