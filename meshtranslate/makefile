CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = parse.o meshtranslate.o 

meshtranslate: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o meshtranslate

lex.c: lex.l
	flex lex.l
	mv lex.yy.c lex.c

parse.c: lex.c parse.y meshtranslate.h
	bison -o parse.c parse.y

parse.o: parse.c meshtranslate.h
	$(CC) $(FLAGS) $(DEFS) -c parse.c

meshtranslate.o: meshtranslate.c meshtranslate.h
	$(CC) $(FLAGS) $(DEFS) -c meshtranslate.c

clean:
	rm -f $(OBJS) lex.c parse.c meshtranslate

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp meshtranslate $(INSTALL_DIR); \
	  cp meshtranslate $(INSTALL_DIR); \
	fi

tar:
	tar -cvf meshtranslate.tar lex.l parse.y meshtranslate.c meshtranslate.h makefile

gzip:
	tar -cvf - lex.l parse.y meshtranslate.c meshtranslate.h makefile | gzip -9c meshtranslate.tar.gz
