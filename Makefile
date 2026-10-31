.PHONY: all audit data validate figures test reproduce clean

VENV := .venv
PYTHON := $(VENV)/bin/python
PYTEST := $(VENV)/bin/pytest

all: reproduce

audit:
	$(PYTHON) -c "import pandas as pd; df=pd.read_csv('data_processed/tcga_prad_recurrence_cohort.csv'); print('Patients:', len(df), 'Events:', df['outcome_recurrence_or_progression'].sum())"

data:
	PYTHONPATH=. $(PYTHON) src/build_recurrence_endpoint.py
	PYTHONPATH=. $(PYTHON) src/cbioportal_endpoint_crosscheck.py

validate:
	PYTHONPATH=. $(PYTHON) src/run_rigorous_validation.py

figures:
	PYTHONPATH=. $(PYTHON) src/run_rigorous_validation.py

test:
	PYTHONPATH=. $(PYTEST) tests/ -v

reproduce: data validate test
	@echo "TCGA-PRAD full scientific pipeline reproduced successfully."

clean:
	rm -rf __pycache__ .pytest_cache
