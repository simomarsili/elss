# change NAME value for a different .pdf filename
NAME=user_guide
$(NAME).pdf: ms.tex ms.bib
#
#default:
	pdflatex ms
	bibtex ms
	pdflatex ms
	pdflatex ms
ifneq ($(NAME), ms)
	mv ms.pdf $(NAME).pdf
endif
	@printf "\ncreated $(NAME).pdf.\n\n"
clean:
	rm -rf ms.log ms.blg ms.bbl ms.aux ?*[~#] 
