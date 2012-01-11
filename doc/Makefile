all: report.pdf

report.pdf: report.tex cg.bib
	rubber --pdf report.tex

clean:
	rubber --clean report.tex
	rm -f report.pdf

.PHONY: clean
