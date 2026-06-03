PYTHON  := python
PACKAGE := puffsat_sim
TESTS   := tests

.PHONY: all run test mypy lint format data clean

all: mypy lint test

run: orekit-data.zip
	$(PYTHON) -m $(PACKAGE).hello_orekit

test:
	pytest $(TESTS)

mypy:
	mypy $(PACKAGE)

lint:
	ruff check $(PACKAGE) $(TESTS)

format:
	ruff format $(PACKAGE) $(TESTS)

# Download the Orekit Earth/time data file (once, ~37 MB).
# Requires gitlab.orekit.org to be reachable.
data: orekit-data.zip

orekit-data.zip:
	$(PYTHON) -c "import orekit_jpype; orekit_jpype.initVM(); \
	    from orekit_jpype.pyhelpers import download_orekit_data_curdir; \
	    download_orekit_data_curdir()"

clean:
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type d -name .mypy_cache -exec rm -rf {} +
	find . -type d -name .ruff_cache -exec rm -rf {} +
	find . -type d -name .pytest_cache -exec rm -rf {} +
	find . -name "*.pyc" -delete
