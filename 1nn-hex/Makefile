.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic .wfc

.qpic.idl: ; $(calcqpi) $*.qpic
.wfc.idl: ; $(mkwavefunctions) $*.wfc

PROJ1=1nnhexfigure.pdf
IDL1=1nn-hex.idl hexwf.idl spf.idl
SRC1=1nn-hex.qpic hexwf.wfc spf.qpic

all: $(PROJ1)

1nn-hex.idl: hexwf.idl 1nn-hex
spf.idl: 1nn-hex

$(PROJ1): $(IDL1)
	python3 MkGraph-1nnhex.py

clean:
	-@rm *.idl *.pdf *.log 2>/dev/null || true
