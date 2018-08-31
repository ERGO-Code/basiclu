# BASICLU/doc/Makefile
#
# The header files of basiclu_ functions contain documentation as C comments.
# For each header file generate a tex file. The tex file embeds the header
# file's content into \begin{verbatim} and \end{verbatim} and removes /* and */.
# The tex files are inputed from userguide.tex.

DOC_FILES=$(patsubst ../include/%.h, %.tex, $(wildcard ../include/basiclu_*.h))

userguide.pdf: $(DOC_FILES) userguide.tex
	pdflatex userguide.tex
	pdflatex userguide.tex

%.tex: ../include/%.h
	/bin/echo '\begin{verbatim}' > $@
	grep -v "/\*\|\*/" $< >> $@
	/bin/echo '\end{verbatim}' >> $@

.PHONY: clean
clean:
	rm -f $(DOC_FILES) userguide.aux userguide.log userguide.toc

.PHONY: purge
purge: clean
	rm -f userguide.pdf
