serial:
	$(MAKE) -C ./src serial

debug:
	$(MAKE) -C ./src debug

hlrn:
	$(MAKE) -C ./src hlrn

test:
	$(MAKE) -C ./test clean

clean:
	$(MAKE) -C ./src clean

testclean:
	$(MAKE) -C ./test	clean

