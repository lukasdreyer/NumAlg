all: compile
PROGNAME1 = ex1
PROGNAME2 = ex2
CFLAGS = -I${PETSC_DIR}/include
SOURCESC = ${PROGNAME1}.c transformation.c utilities.c griddata.c  #$(wildcard *.c)
OBJ = $(SOURCESC:.c=.o)
CLEANFILES = ${OBJ} ${PROGNAME1}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

compile:${OBJ}
	-${CLINKER} -o ${PROGNAME1} ${OBJ} ${PETSC_SYS_LIB}


runex1: 
	-@${MPIEXEC} -n 1 ./${PROGNAME1}

help:
	./${PROGNAME1} -help intro
	
clean ::
	-${RM} ${CLEANFILES}

.PHONY: all compile runtex1 help clean

