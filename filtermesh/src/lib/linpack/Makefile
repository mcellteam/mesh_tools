# Converted fortran code to c using "f2c -c"
# Found file r_sign.c in f2c library and created f2chelper.c with it.
# On alpha, had to modify 'long' to be 'int' in f2c.h

CFLAGS = -O3 -march=pentiumpro

# FTN = f77
# FTNFLAGS = -O
# .SUFFIXES: .f .o
# .f.o:
# 	$(FTN) $(FTNFLAGS) -c $<

LIB = liblinpack.a
OBJ = saxpy.o scopy.o sdot.o snrm2.o sqrdc.o sqrsl.o \
	sscal.o sswap.o f2chelper.o
all:	${LIB}
test:	testlinpack
clean:
	rm -f *.o *.a
${LIB}:	${OBJ}
	ar r ${LIB} ${OBJ}
	ranlib ${LIB}
# -L for sparc
# testlinpack: testlinpack.o ${LIB}
# 	$(CC) $(CFLAGS) $(LDFLAGS) $@.o ${LIB} -o $@ \
# 	  -L/usr/lang/SC1.0.1 -lF77 -lm
testlinpack: testlinpack.o ${LIB}
	$(CC) $(CFLAGS) $(LDFLAGS) $@.o ${LIB} -o $@ -lm
