.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic

.qpic.idl: ; $(calcqpi) $*.qpic

PROJ1=1nn-dx2figure.pdf
IDL1=1nn-dx2j.idl 1nn-dx2jn.idl
SRC1=1nn-dx2j.qpic 1nn-dx2jn.qpic

1nn-dx2j.idl: 1nn-dx2
1nn-dx2jn.idl: 1nn-dx2

all: $(PROJ1)

$(PROJ1): $(IDL1)
	python3 MkPlot-josephson.py

clean:
	-@rm *.idl *.pdf *.log 2>/dev/null || true
