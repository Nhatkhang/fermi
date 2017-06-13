#
# makefile to compile "control"
#
#  Authors: 
#
#  Miguel Zavala
#  Guido Giuntoli
#

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
	mpicxx ${CFLAGS} -o $@ $< ${LIBS} -lm 

clean_:	    
	rm -f control
