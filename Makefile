.DEFAULT_GOAL := all

all:
	cd 1nn; $(MAKE) all
	cd 1nnz; $(MAKE) all
	cd 1nn-sc; $(MAKE) all
	cd 1nn-josephson; $(MAKE) all
	cd 1nn-hex; $(MAKE) all
	cd rashba; $(MAKE) all
	cd topo; $(MAKE) all
	cd ssh; $(MAKE) all

default:
	$(MAKE) all

install:
	cd 1nn; $(MAKE) install
	cd 1nnz; $(MAKE) install
	cd 1nn-sc; $(MAKE) install
	cd 1nn-josephson; $(MAKE) install
	cd 1nn-hex; $(MAKE) install
	cd rashba; $(MAKE) install
	cd topo; $(MAKE) install
	cd ssh; $(MAKE) install

depend:
	cd 1nn; $(MAKE) depend
	cd 1nnz; $(MAKE) depend
	cd 1nn-sc; $(MAKE) depend
	cd 1nn-josephson; $(MAKE) depend
	cd 1nn-hex; $(MAKE) depend
	cd rashba; $(MAKE) depend
	cd topo; $(MAKE) depend
	cd ssh; $(MAKE) depend

clean:
	cd 1nn; $(MAKE) clean
	cd 1nnz; $(MAKE) clean
	cd 1nn-sc; $(MAKE) clean
	cd 1nn-josephson; $(MAKE) clean
	cd 1nn-hex; $(MAKE) clean
	cd rashba; $(MAKE) clean
	cd topo; $(MAKE) clean
	cd ssh; $(MAKE) clean
