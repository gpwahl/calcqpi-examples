.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic .wfc

.qpic.idl: ; $(calcqpi) $*.qpic
.wfc.idl: ; $(mkwavefunctions) $*.wfc

PROJ1=1nnhexfigure.pdf
IDL1=1nn-hex.idl hexwf.idl bshist-hex.idl
SRC1=1nn-hex.qpic hexwf.wfc

all: $(PROJ1)

1nn-hex.idl: hexwf.idl

bshist-hex.idl:
	writebandstructure -obshist-hex.idl -l41 -n200 -r-0.03,0.03 1nn-hex

$(PROJ1): $(IDL1)
	python3 MkGraph-1nnhex.py

clean:
	-rm *.idl *.pdf *.log
