# The main makefile
#

all: src

clean: force-look
	@echo ==== Cleaning package ====
	(cd src; make clean)

src: force-look
	(cd src; make all)

force-look:
	@true

cleanmats:
	rm *.nb
	rm *.emm
	rm *.log
