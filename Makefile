serial:
	$(MAKE) -C ./src serial

debug:
	$(MAKE) -C ./src debug

mpi:
	$(MAKE) -C ./src mpi

.PHONY: test
test:
	$(MAKE) -C ./test test

clean:
	$(MAKE) -C ./src clean

testclean:
	$(MAKE) -C ./test	clean

PREFIX ?= /usr/local

install:
	mkdir -p $(PREFIX)/bin
	# serial programs
	cp bin/rixs-pathway-serial $(PREFIX)/bin/ || true
	cp bin/rixs-oscstr-serial $(PREFIX)/bin/ || true
	cp bin/rixs-coherence-serial $(PREFIX)/bin/ || true
