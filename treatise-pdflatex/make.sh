#!/usr/bin/env bash

NAME=main

clear
#rm ${NAME}.out ${NAME}.log ${NAME}.bcf ${NAME}.aux ${NAME}.run.xml
pdflatex -halt-on-error ${NAME}.tex
biber ${NAME}
pdflatex -halt-on-error ${NAME}.tex
cp -f ${NAME}.pdf ${NAME}_1.pdf
cp -f ${NAME}.pdf ${NAME}_2.pdf


