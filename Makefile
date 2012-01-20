# The main makefile
#

all: src doc

doc: force-look
	@echo ==== Making documentation ====
	(cd doc; make)

clean: force-look
	@echo ==== Cleaning package ====
	(cd src; make clean)
	(cd doc; make clean)

src: force-look
	@echo ==== Building tools ====
	(cd src; make all)

force-look:
	@true

cleanmats:
	rm mat-check-*.nb
	rm *.emm
	rm *.log
