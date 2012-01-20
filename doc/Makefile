all: report.pdf

report.pdf: report.tex cg.bib img/*
	rubber --pdf report.tex

clean:
	rubber --clean report.tex
	rm -vf report.pdf

.PHONY: clean
