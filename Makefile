prefix := /SNS/REF_L/shared

all:
	@echo "Run 'make install' to install the automated reduction code for LR"


install:
	cp -R scripts/autoreduce/*.py $(prefix)/autoreduce
	cp -R scripts/shared/*.py $(prefix)


.PHONY: install
