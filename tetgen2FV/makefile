CC = gcc
DEFS =
#CFLAGS = -Wall -O3 
#CFLAGS = -Wall -O3
CFLAGS = -Wall -ggdb3 -O0 -pg
LIBS = -lm
# LIBS = -lefence -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = tetgennodeparse.o tetgenfaceparse.o tetgenelemparse.o output_dx.o output_edsim.o output_renderman.o strfunc.o hash.o c_3d.o polyarea.o polyvol.o polycent.o vector.o voronoi.o tetgen2FV.o 

tetgen2FV: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o tetgen2FV

tetgennodelex.c: tetgennodelex.l
	flex -Cfra -Ptetgennode -otetgennodelex.c tetgennodelex.l

tetgennodeparse.c: tetgennodelex.c tetgennodeparse.y tetgen2FV.h
	bison -p tetgennode -o tetgennodeparse.c tetgennodeparse.y

tetgennodeparse.o: tetgennodeparse.c tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c tetgennodeparse.c

tetgenfacelex.c: tetgenfacelex.l
	flex -Cfra -Ptetgenface -otetgenfacelex.c tetgenfacelex.l

tetgenfaceparse.c: tetgenfacelex.c tetgenfaceparse.y tetgen2FV.h
	bison -p tetgenface -o tetgenfaceparse.c tetgenfaceparse.y

tetgenfaceparse.o: tetgenfaceparse.c tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c tetgenfaceparse.c

tetgenelemlex.c: tetgenelemlex.l
	flex -Cfra -Ptetgenelem -otetgenelemlex.c tetgenelemlex.l

tetgenelemparse.c: tetgenelemlex.c tetgenelemparse.y tetgen2FV.h
	bison -p tetgenelem -o tetgenelemparse.c tetgenelemparse.y

tetgenelemparse.o: tetgenelemparse.c tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c tetgenelemparse.c

output_dx.o: output_dx.c output_dx.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c output_dx.c

output_edsim.o: output_edsim.c output_edsim.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c output_edsim.c

output_renderman.o: output_renderman.c output_renderman.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c output_renderman.c

strfunc.o: strfunc.c strfunc.h
	$(CC) $(CFLAGS) $(DEFS) -c strfunc.c

hash.o: hash.c hash.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c hash.c

c_3d.o: c_3d.c c_3d.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c c_3d.c

polyarea.o: polyarea.c polyarea.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c polyarea.c

polyvol.o: polyvol.c polyvol.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c polyvol.c

polycent.o: polycent.c polycent.h vector.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c polycent.c

vector.o: vector.c vector.h
	$(CC) $(CFLAGS) $(DEFS) -c vector.c

voronoi.o: voronoi.c voronoi.h c_3d.h tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c voronoi.c

tetgen2FV.o: tetgen2FV.c tetgen2FV.h
	$(CC) $(CFLAGS) $(DEFS) -c tetgen2FV.c

clean:
	rm -f $(OBJS) tetgennodelex.c tetgennodeparse.c tetgenfacelex.c tetgenfaceparse.c tetgenelemlex.c tetgenelemparse.c tetgen2FV

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp tetgen2FV $(INSTALL_DIR); \
	  cp tetgen2FV $(INSTALL_DIR); \
	fi

tar:
	tar -cvf tetgen2FV.tar lex.l parse.y tetgen2FV.c tetgen2FV.h makefile

gzip:
	tar -cvf - lex.l parse.y tetgen2FV.c tetgen2FV.h makefile | gzip -9c tetgen2FV.tar.gz
