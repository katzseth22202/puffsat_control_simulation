PYTHON  := python
PACKAGE := puffsat_sim
TESTS   := tests

.PHONY: all run capstone test test-integration mypy lint format format-check data clean

all: mypy lint format-check test

run: orekit-data.zip
	$(PYTHON) -m $(PACKAGE).truth_model

# Open-loop dispersion capstone — smoke-sized ensemble (N=50; the resolved-tail
# result needs N=10^3-10^4, a longer job).
capstone: orekit-data.zip
	$(PYTHON) -m $(PACKAGE).montecarlo

# Pure-Python unit tests — no JVM required.
test:
	pytest $(TESTS) --ignore=$(TESTS)/integration

# Integration tests — require a live Orekit JVM and orekit-data.zip in the cwd.
test-integration: orekit-data.zip
	pytest -m integration $(TESTS)/integration

mypy:
	mypy $(PACKAGE)

lint:
	ruff check $(PACKAGE) $(TESTS)

format:
	ruff format $(PACKAGE) $(TESTS)

# Gate: fail if any file is not ruff-formatted (run by `make all`).
format-check:
	ruff format --check $(PACKAGE) $(TESTS)

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
