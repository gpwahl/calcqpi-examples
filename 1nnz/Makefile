.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic

.qpic.idl: ; $(calcqpi) $*.qpic

PROJ1=1nnzfigure.pdf
IDL1=nn1zqpi.idl spf-bulk.idl
SRC1=nn1zqpi.qpic spf-bulk.qpic

nn1zqpi.idl: 1nn-z
spf-bulk.idl: 1nn-z

all: $(PROJ1)

$(PROJ1): $(IDL1)
	python3 MkPlot-1nnz.py

clean:
	-@rm *.idl *.pdf *.log 2>/dev/null || true
