all: compile1 compile2

PROGNAME1 = ex1
PROGNAME2 = ex2
CFLAGS = -I${PETSC_DIR}/include
SOURCESADD = transformation.c utilities.c griddata.c matshell.c  #$(wildcard *.c)
SOURCESP1 = ${PROGNAME1}.c ${SOURCESADD}
SOURCESP2 = ${PROGNAME2}.c ${SOURCESADD}
SOURCES =  ${PROGNAME1}.c ${PROGNAME2}.c ${SOURCESADD}	
OBJP1 = $(SOURCESP1:.c=.o)
OBJP2 = $(SOURCESP2:.c=.o)
OBJS = $(SOURCES:.c=.o)
CLEANFILES = ${OBJ} ${PROGNAME1} ${PROGNAME2}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

compile1:${OBJP1}
	-${CLINKER} -o ${PROGNAME1} ${OBJP1} ${PETSC_SYS_LIB}

compile2:${OBJP2}
	-${CLINKER} -o ${PROGNAME2} ${OBJP2} ${PETSC_SYS_LIB}	

runex1: 
	-@${MPIEXEC} -n 1 ./${PROGNAME1}

runex2: 
	-@${MPIEXEC} -n 1 ./${PROGNAME2}

	
help:
	./${PROGNAME1} -help intro
	
clean ::
	-${RM} ${CLEANFILES}

.PHONY: all compile1 compile2 runex1 runex2 help clean

