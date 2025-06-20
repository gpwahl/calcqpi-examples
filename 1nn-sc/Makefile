.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic

.qpic.idl: ; $(calcqpi) $*.qpic

PROJ1=1nn-dx2figure.pdf
IDL1=1nn.idl 1nn-dx2.idl
SRC1=1nn.qpic 1nn-dx2.qpic

all: $(PROJ1)

$(PROJ1): $(IDL1)
	python3 MkPlot-1nn-dx2.py

clean:
	-rm *.idl *.pdf *.log
