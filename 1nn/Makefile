.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic

.qpic.idl: ; $(calcqpi) $*.qpic

PROJ1=1nnfigure.pdf
IDL1=1nnmodel.idl spf.idl
SRC1=1nnmodel.qpic spf.qpic

all: $(PROJ1)

$(PROJ1): $(IDL1)
	python3 MkPlot-1nn.py

clean:
	-@rm *.idl *.pdf *.log 2>/dev/null || true
