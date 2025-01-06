# bash sell to correctly interpret the double brackets in the conditions below
SHELL=/bin/bash
# https://www.gnu.org/software/make/manual/html_node/One-Shell.html
# Required to prevent having to use "conda init"

# all the lines in a recipe are passed to a single invocation of the shell.
.ONESHELL:

# list of all phony targets, alphabetically sorted
.PHONY: help conda docs test install

help:
    # This nifty perl one-liner collects all commnents headed by the double "#" symbols next to each target and recycles them as comments
	@perl -nle'print $& if m{^[a-zA-Z_-]+:.*?## .*$$}' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-25s\033[0m %s\n", $$1, $$2}'


PREFIX := /SNS/REF_L/shared
install:  ## Install the automated reduction code for LR
	versioningit -w
	cp -R scripts/livereduce/*.py $(PREFIX)/livereduce
	cp -R scripts/autoreduce/*.py $(PREFIX)/autoreduce
	cp -R scripts/shared/*.py $(PREFIX)
	cp -R scripts/shared/.*.conf $(PREFIX)
	cp -R launcher $(PREFIX)
	cp -R reduction $(PREFIX)
	cp -R xrr/xrr_processing.py $(PREFIX)

# Note that the extra activate is needed to ensure that the activate floats env to the front of PATH
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate
conda-env:  ## creates conda environment `lr_reduction` and installs package `lr_reduction` in editable mode
	conda env create --solver=libmamba --file ./environment.yml
	$(CONDA_ACTIVATE) lr_reduction
	pip install -e .

docs: ## Build the documentation
	mkdocs build

test-all:  ## Run all tests
	pytest ./test
