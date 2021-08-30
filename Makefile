prefix := /SNS/REF_L/shared

all:
	@echo "Run 'make install' to install the automated reduction code for LR"


install:
	cp -R scripts/autoreduce/*.py $(prefix)/autoreduce
	cp -R scripts/shared/*.py $(prefix)
	cp -R scripts/shared/.*.conf $(prefix)
	cp -R 30Hz/template_reduction.py $(prefix)
	cp -R 30Hz/time-resolved-ui.py $(prefix)
	cp -R xrr/xrr_processing.py $(prefix)

.PHONY: install
