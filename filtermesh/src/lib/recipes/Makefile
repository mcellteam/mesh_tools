BASEDIR = ../..
OPTIM = -O3 -march=pentiumpro
INC = -I${BASEDIR}/include
CFLAGS = ${OPTIM} ${INC}

LIB = librecipes.a
OBJ = \
	nrutil.o ludcmp.o lubksb.o svdcmp.o svbksb.o jacobi.o eigsrt.o
all:	${LIB}
clean:
	rm -f *.o *.a
${LIB}:	${OBJ}
	ar r ${LIB} ${OBJ}
	ranlib ${LIB}
