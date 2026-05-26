TEXFILES = main.tex supplementary.tex cover_letter.tex ../zotero.bib

PDFLATEX = pdflatex -interaction=nonstopmode

all: main.pdf supplementary.pdf cover_letter.pdf

main.pdf: main.tex ../zotero.bib
	$(PDFLATEX) main.tex
	bibtex main
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex

supplementary.pdf: supplementary.tex ../zotero.bib
	$(PDFLATEX) supplementary.tex
	bibtex supplementary
	$(PDFLATEX) supplementary.tex
	$(PDFLATEX) supplementary.tex

cover_letter.pdf: cover_letter.tex
	$(PDFLATEX) cover_letter.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.pdf *.run.xml *.toc

.PHONY: all clean
