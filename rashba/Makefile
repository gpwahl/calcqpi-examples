.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic

.qpic.idl: ; $(calcqpi) $*.qpic

PROJ1=rashba.pdf
IDL1=qpi1nnrashba.idl spf.idl
SRC1=qpi1nnrashba.qpic spf.qpic

all: $(PROJ1)

$(PROJ1): $(IDL1)
	python3 MkGraph-rashba.py

clean:
	-rm *.idl *.pdf *.log
