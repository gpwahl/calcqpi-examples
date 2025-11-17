.DEFAULT_GOAL := all

include ../config.mk

.SUFFIXES: .idl .qpic

.qpic.idl: ; $(calcqpi) $*.qpic

PROJ1=qpissh.pdf
IDL1=qpibulk.idl qpisurface.idl spfbulk.idl spfsurface.idl
SRC1=qpibulk.qpic qpisurface.qpic spfbulk.qpic spfsurface.qpic

qpibulk.idl: ssh1nn.dat
qpisurface.idl: ssh1nn.dat
spfbulk.idl: ssh1nn.dat
spfsurface.idl: ssh1nn.dat

all: $(PROJ1)

$(PROJ1): $(IDL1)
	python3 MkGraph-ssh.py

clean:
	-@rm *.idl *.pdf *.log 2>/dev/null || true
