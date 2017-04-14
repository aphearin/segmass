set local=./

latex segmass.tex
bibtex segmass
latex segmass.tex
latex segmass.tex
dvips segmass.dvi -Ppdf -G0 -z -t a4
open segmass.ps
