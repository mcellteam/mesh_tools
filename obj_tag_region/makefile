CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.bison.o vector.o obj_tag_region.o 

obj_tag_region: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o obj_tag_region

lex.flex.c: lex.l
	flex -Cfra lex.l

parse.bison.c: lex.flex.c parse.y obj_tag_region.h
	bison parse.y

parse.bison.o: parse.bison.c obj_tag_region.h
	$(CC) $(FLAGS) $(DEFS) -c parse.bison.c

vector.o: vector.h
	$(CC) $(FLAGS) $(DEFS) -c vector.c

obj_tag_region.o: obj_tag_region.c obj_tag_region.h
	$(CC) $(FLAGS) $(DEFS) -c obj_tag_region.c

clean:
	rm -f $(OBJS) lex.flex.c parse.bison.c obj_tag_region

install: obj_tag_region
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp obj_tag_region $(INSTALL_DIR); \
	  cp obj_tag_region $(INSTALL_DIR); \
	fi

tar:
	tar -cvf obj_tag_region.tar lex.l parse.y obj_tag_region.c obj_tag_region.h makefile

gzip:
	tar -cvf - lex.l parse.y obj_tag_region.c obj_tag_region.h makefile | gzip -9c obj_tag_region.tar.gz
