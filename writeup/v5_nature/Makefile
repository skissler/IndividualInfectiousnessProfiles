PDFLATEX = pdflatex -interaction=nonstopmode

all: main.pdf supplementary.pdf cover_letter.pdf

main.pdf: main.tex ../zotero.bib
	-$(PDFLATEX) main.tex
	-bibtex main
	-$(PDFLATEX) main.tex
	-$(PDFLATEX) main.tex

science: main-science.pdf

main-science.pdf: main.tex science-style.sty ../zotero.bib
	-$(PDFLATEX) -jobname=main-science '\def\sciencemode{}\input{main}'
	-bibtex main-science
	-$(PDFLATEX) -jobname=main-science '\def\sciencemode{}\input{main}'
	-$(PDFLATEX) -jobname=main-science '\def\sciencemode{}\input{main}'

supplementary.pdf: supplementary.tex ../zotero.bib
	-$(PDFLATEX) supplementary.tex
	-bibtex supplementary
	-$(PDFLATEX) supplementary.tex
	-$(PDFLATEX) supplementary.tex

cover_letter.pdf: cover_letter.tex
	-$(PDFLATEX) cover_letter.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.pdf *.run.xml *.toc

.PHONY: all science clean
