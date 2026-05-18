TEXFILES = main.tex authors.tex abstract.tex introduction.tex model.tex uncontrolled.tex \
           controlled.tex identifiability.tex discussion.tex acknowledgements.tex \
           supplementary.tex zotero.bib

PDFLATEX = pdflatex -interaction=nonstopmode

all: main.pdf

main.pdf: $(TEXFILES)
	$(PDFLATEX) main.tex
	bibtex main
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex

dark: $(TEXFILES)
	$(PDFLATEX) "\def\darkmode{}\input{main}" -jobname=main
	bibtex main
	$(PDFLATEX) "\def\darkmode{}\input{main}" -jobname=main
	$(PDFLATEX) "\def\darkmode{}\input{main}" -jobname=main

clean:
	rm -f main.aux main.bbl main.blg main.log main.out main.pdf main.run.xml

.PHONY: all dark clean
