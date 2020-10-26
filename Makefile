serial:
	$(MAKE) -C ./src serial

debug:
	$(MAKE) -C ./src debug

hlrn:
	$(MAKE) -C ./src hlrn

.PHONY: test
test:
	$(MAKE) -C ./test test

clean:
	$(MAKE) -C ./src clean

testclean:
	$(MAKE) -C ./test	clean

