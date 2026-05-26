TEXFILES = main.tex ../zotero.bib

PDFLATEX = pdflatex -interaction=nonstopmode

all: main.pdf

main.pdf: $(TEXFILES)
	$(PDFLATEX) main.tex
	bibtex main
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.pdf *.run.xml

.PHONY: all clean
