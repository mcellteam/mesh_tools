CC = gcc
DEFS =
FLAGS = -O3 -march=pentiumpro
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o vector.o mesh_tag_region.o 

mesh_tag_region: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o mesh_tag_region

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y mesh_tag_region.h
	bison -o parse.c parse.y

parse.o: parse.c mesh_tag_region.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

vector.o: vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

mesh_tag_region.o: mesh_tag_region.c mesh_tag_region.h
	$(CC) $(FLAGS) $(DEFS) -c mesh_tag_region.c

clean:
	rm -f $(OBJS) lex.c parse.c mesh_tag_region

install: mesh_tag_region
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp mesh_tag_region $(INSTALL_DIR); \
	  cp mesh_tag_region $(INSTALL_DIR); \
	fi

tar:
	tar -cvf mesh_tag_region.tar lex.l parse.y mesh_tag_region.c mesh_tag_region.h makefile

gzip:
	tar -cvf - lex.l parse.y mesh_tag_region.c mesh_tag_region.h makefile | gzip -9c mesh_tag_region.tar.gz
