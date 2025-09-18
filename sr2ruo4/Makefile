.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic

.qpic.idl: ; $(calcqpi) $*.qpic

PROJ1=214figure.pdf
IDL1=214qpi01.idl 214qpi10.idl
SRC1=214qpi01.qpic 214qpi10.qpic

all: $(PROJ1)

$(PROJ1): $(IDL1)
	python3 MkPlot-214.py

clean:
	-@rm *.idl *.pdf *.log 2>/dev/null || true
