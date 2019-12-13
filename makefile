all: compile1 compile2 compileut

PROGNAME1 = ex1
PROGNAME2 = ex2
UNITTEST = unittest
CFLAGS = -I${PETSC_DIR}/include
SOURCESADD = transformation.c utilities.c griddata.c matshell.c mat2x2.c testfunctions.c solver.c #$(wildcard *.c)
SOURCESP1 = ${PROGNAME1}.c ${SOURCESADD}
SOURCESP2 = ${PROGNAME2}.c ${SOURCESADD}
SOURCESUT = ${UNITTEST}.c ${SOURCESADD}
SOURCES =  ${PROGNAME1}.c ${PROGNAME2}.c ${UNITTEST}.c ${SOURCESADD}	
OBJP1 = $(SOURCESP1:.c=.o)
OBJP2 = $(SOURCESP2:.c=.o)
OBJUT = $(SOURCESUT:.c=.o)
OBJS = $(SOURCES:.c=.o)
CLEANFILES = ${OBJ} ${PROGNAME1} ${PROGNAME2} ${UNITTEST}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

compile1:${OBJP1}
	-${CLINKER} -o ${PROGNAME1} ${OBJP1} ${PETSC_SYS_LIB}

compile2:${OBJP2}
	-${CLINKER} -o ${PROGNAME2} ${OBJP2} ${PETSC_SYS_LIB}	

compileut:${OBJUT}
	-${CLINKER} -g -o ${UNITTEST} ${OBJUT} ${PETSC_SYS_LIB}

debug2:${OBJP2}
	-${CLINKER} -g -o ${PROGNAME2} ${OBJP2} ${PETSC_SYS_LIB}	

runex1: 
	-@${MPIEXEC} -n 1 ./${PROGNAME1}

runex2: 
	-@${MPIEXEC} -n 1 ./${PROGNAME2}
	
unittest: 
	-@${MPIEXEC} -n 1 ./${UNITTEST}
	
help:
	./${PROGNAME1} -help intro
	
clean ::
	-${RM} ${CLEANFILES}

.PHONY: all compile1 compile2 compileut runex1 runex2 debug 2 unittest help clean

