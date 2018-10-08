LX = xelatex
MI = makeindex
BX = biber

FILEN = article
TEXFLAGS =
PREFLAGS = -no-pdf

once:
	-$(LX) $(TEXFLAGS) $(FILEN) 

quick:
	-$(LX) $(TEXFLAGS) $(PREFLAGS) $(FILEN)
	printf "\033c"
	@echo
	@echo " === LAST RUN === "
	@echo
	-$(LX) $(TEXFLAGS) $(FILEN) | grep 'Warning\|Error'

clean:
	\rm *.adx *.and *.aux *.bbl *.blg *.idx *.ilg *.ldx *.lnd *.log *.out *.pdf *.rdx *.run.xml *.sdx *.snd *.toc *.wdx *.xdv
	\rm chapters/*.aux

all:
	-$(LX) $(TEXFLAGS) $(PREFLAGS) $(FILEN)
	-$(BX) $(FILEN)
	-$(LX) $(TEXFLAGS) $(PREFLAGS) $(FILEN)
	printf "\033c"
	@echo
	@echo " === LAST RUN === "
	@echo
	-$(LX) $(TEXFLAGS) $(FILEN) | grep 'Warning\|Error'


view:
	/Applications/Skim.app/Contents/MacOS/Skim $(FILEN).pdf & 

edit:
	mvim -c ':set spell spelllang=en' -c ':nnoremap <F15> ]s' -c ':nnoremap <F14> [s' article.tex
