CC = gcc
DEFS =
FLAGS = -O3
LIBS = -lm
INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh
OBJS = netgenparse.bison.o strfunc.o netgen2obj.o 

netgen2obj: $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o netgen2obj

netgenlex.flex.c: netgenlex.l
	flex -Cfra -P netgen netgenlex.l

netgenparse.bison.c: netgenlex.flex.c netgenparse.y netgen2obj.h
	bison netgenparse.y

netgenparse.bison.o: netgenparse.bison.c netgen2obj.h
	$(CC) $(FLAGS) $(DEFS) -c netgenparse.bison.c

strfunc.o: strfunc.c strfunc.h
	$(CC) $(FLAGS) $(DEFS) -c strfunc.c

netgen2obj.o: netgen2obj.c netgen2obj.h
	$(CC) $(FLAGS) $(DEFS) -c netgen2obj.c

clean:
	rm -f $(OBJS) netgenlex.flex.c netgenparse.bison.c netgen2obj

install:
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp netgen2obj $(INSTALL_DIR); \
	  cp netgen2obj $(INSTALL_DIR); \
	fi

tar:
	tar -cvf netgen2obj.tar netgenlex.l netgenparse.y netgen2obj.c netgen2obj.h strfunc.c strfunc.h makefile

gzip:
	tar -cvf - netgenlex.l netgenparse.y netgen2obj.c netgen2obj.h strfunc.c strfunc.h makefile | gzip -9c netgen2obj.tar.gz
