.DEFAULT_GOAL := all

SUBDIRS := 1nn 1nnz 1nn-sc 1nn-josephson 1nn-hex rashba topo ssh sr2ruo4

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ 

default:
	$(MAKE) all

clean:
	-@for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done

.PHONY: all $(SUBDIRS)
