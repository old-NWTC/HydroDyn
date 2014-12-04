# 2014.03.28 -- ADP
#  This makefile is for compiling the guide for the WAMIT2 development documentation (2nd order hydrodyn)



Mainfile:=\
   HD-2ndOrd--implimentation.tex

Result:=$(Mainfile:.tex=.pdf)


Secondary:= \
	chaps/Chap.Second_order_calcutions.tex						\
	chaps/Chap.Background.tex										\
	chaps/Chap.Bi-Directional_Waves.tex							\
	chaps/Chap.Development_schedule.tex							\
	chaps/Chap.Interpolation_Algorithms.tex					\
	chaps/Chap.Overview.tex											\
	chaps/Chap.Second_order_force_calculations.tex			\
	chaps/Chap.WAMIT2_module_Architecture.tex					\
	chaps/Chap.WAMIT2_module_Integration.tex

figures:=   \
	./chaps/figures/HydroDyn_2nd_order--WAMIT2_Init.pdf   \
	./chaps/figures/WAMIT2_Dev.pdf                        \
	./chaps/figures/HydroDyn_2nd_order--WAMIT2_Init.pdf   \
	./chaps/figures/WAMIT2_Dev.pdf


default: 
	pdflatex -shell-escape $(Mainfile)

all:
	pdflatex -shell-escape $(Mainfile)
	bibtex $(Mainfile:.tex=.aux)
	pdflatex -shell-escape $(Mainfile)
	pdflatex -shell-escape $(Mainfile)

clean:
	-rm *~
	-rm */*~
	-rm */*/*~

cleanall: clean
	-rm *aux *lof *log *lot *out *pyg *toc *bbl *question *blg *tdo
	#-rm $(Mainfile:.tex=.pdf)
