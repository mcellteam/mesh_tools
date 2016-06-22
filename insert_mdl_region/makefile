INSTALL_DIR = $(HOME)/bin

SHELL = /bin/sh

install: insert_mdl_region.py
	@if [ "$(INSTALL_DIR)" ]; then \
	  if [ ! -s $(INSTALL_DIR) ]; then \
	    echo mkdir $(INSTALL_DIR); \
	    mkdir $(INSTALL_DIR); \
	  fi; \
	  echo cp insert_mdl_region.py $(INSTALL_DIR); \
	  cp insert_mdl_region.py $(INSTALL_DIR); \
	fi

