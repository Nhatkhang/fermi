#
# makefile to compile "control"
#
#  Authors: 
#
#  Miguel Zavala
#  Guido Giuntoli
#

PWD:= $(shell pwd)

CFLAGS=-g -O0

ROOT     = ../..
LIBPLEPP = $(ROOT)/lib/PLEPP/LIBPLEPP
LIBPLE   = $(ROOT)/lib/PLEPP/LIBPLE
CFLAGS  += -DCOMMDOM=$(COMMDOM)
CFLAGS  += -I$(LIBPLEPP)/Include/
CFLAGS  += -I$(LIBPLEPP)/Wrappers/CC/
LIBS     = $(LIBPLEPP)/Wrappers/CC/libcommdom.a 
LIBS    += $(LIBPLE)/Execs/lib/libple.a 

.PHONY: clean_

all: control 
	
control: control.c
	${PETSC_COMPILE} ${CFLAGS} -o $@ $< 

clean_:	    
	rm -f control

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
