.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic

.qpic.idl: ; $(calcqpi) $*.qpic

PROJ1=qpissh.pdf
IDL1=qpibulk.idl qpisurface.idl spfbulk.idl spfsurface.idl
SRC1=qpibulk.qpic qpisurface.qpic spfbulk.qpic spfsurface.qpic

qpibulk.idl: topo.dat
qpisurface.idl: topo.dat
spfbulk.idl: topo.dat
spfsurface.idl: topo.dat

all: $(PROJ1)

$(PROJ1): $(IDL1)
	python3 MkGraph-topo.py

clean:
	-@rm *.idl *.pdf *.log 2>/dev/null || true
