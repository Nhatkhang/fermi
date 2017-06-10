#############################################################
#
# FERMI makefile
#
# Authors: 
# Miguel Zavala
# Guido Giuntoli
#

PWD:= $(shell pwd)

SRC_DIR=./src
OBJ_DIR=./obj
DEP_DIR=${PWD}/inc

CFLAGS=-g -O0
	
DEPS = ${DEP_DIR}/fermi.h              \
       ${DEP_DIR}/gmsh.h               \
       ${DEP_DIR}/mesh.h               \
       ${DEP_DIR}/fem.h                \
       ${DEP_DIR}/fun.h                \
       ${DEP_DIR}/list.h               \
       ${DEP_DIR}/types.h              \
       ${DEP_DIR}/global.h             \
       ${DEP_DIR}/utils.h              

OBJ  = ${OBJ_DIR}/fer_main.o           \
       ${OBJ_DIR}/mesh.o               \
       ${OBJ_DIR}/list.o               \
       ${OBJ_DIR}/gmsh.o               \
       ${OBJ_DIR}/fem.o                \
       ${OBJ_DIR}/fun.o                \
       ${OBJ_DIR}/fer_assm.o           \
       ${OBJ_DIR}/fer_init.o           \
       ${OBJ_DIR}/fer_finish.o         \
       ${OBJ_DIR}/fer_boun.o           \
       ${OBJ_DIR}/fer_step.o           \
       ${OBJ_DIR}/fer_elem.o           \
       ${OBJ_DIR}/fer_solv.o           \
       ${OBJ_DIR}/fer_power.o          \
       ${OBJ_DIR}/fer_rods.o           \
       ${OBJ_DIR}/lst2msh.o            \
       ${OBJ_DIR}/utils.o              \
       ${OBJ_DIR}/fer_output.o         \
       ${OBJ_DIR}/fer_parser.o             

.PHONY: clean_

PLEPP=0

ifneq ($(PLEPP), 0)

  ROOT      = $(shell pwd)

  LIBPLEPP  = $(ROOT)/lib/PLEPP/LIBPLEPP
  LIBPLE    = $(ROOT)/lib/PLEPP/LIBPLE

  CFLAGS   += -DCOMMDOM=$(COMMDOM)
  CFLAGS   += -I$(LIBPLEPP)/Include
  CFLAGS   += -I$(LIBPLEPP)/Wrappers/CC

  LIBS      = $(LIBPLEPP)/Wrappers/CC/libcommdom.a 
  LIBS     += $(LIBPLE)/Execs/lib/libple.a 

endif

all: ${OBJ} 
	mpicc -o fermi  $^ $(LIBS) ${SLEPC_EPS_LIB} 
	
${OBJ_DIR}/%.o: ${SRC_DIR}/%.c $(DEPS) 
	${PETSC_COMPILE} -c ${CFLAGS} -o $@ $< -I${DEP_DIR}

clean_:	    
	rm -f $(OBJ) fermi

vars:
	@echo $(LIBS)

include ${PETSC_DIR}/lib/petsc/conf/variables	
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
