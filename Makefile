prefix := /SNS/REF_L/shared

all:
	@echo "Run 'make install' to install the automated reduction code for LR"


install:
	cp -R scripts/autoreduce/*.py $(prefix)/autoreduce
	cp -R scripts/shared/*.py $(prefix)
	cp -R scripts/shared/.*.conf $(prefix)
	cp -R launcher $(prefix)
	cp -R reduction $(prefix)
	cp -R xrr/xrr_processing.py $(prefix)

.PHONY: install
