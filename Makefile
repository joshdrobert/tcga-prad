.PHONY: all setup audit data validate figures test smoke reproduce oral-sprint clean

VENV := .venv
PYTHON := $(VENV)/bin/python
PYTEST := $(VENV)/bin/pytest

all: reproduce

setup:
	@$(PYTHON) -c "import numpy, scipy, pandas, sklearn, lifelines; print('Environment verified.')"

audit:
	$(PYTHON) scripts/audit_conference_claims.py

oral-sprint:
	PYTHONPATH=. $(PYTHON) src/run_oral_sprint_analysis.py
	$(PYTHON) scripts/audit_oral_candidate_claims.py
	PYTHONPATH=. $(PYTEST) tests/test_oral_candidate_abstract.py -v
	@echo "Oral candidate sprint pipeline verified successfully."

data:
	PYTHONPATH=. $(PYTHON) src/build_recurrence_endpoint.py
	PYTHONPATH=. $(PYTHON) src/cbioportal_endpoint_crosscheck.py

validate:
	PYTHONPATH=. $(PYTHON) src/run_rigorous_validation.py

figures:
	PYTHONPATH=. $(PYTHON) src/run_rigorous_validation.py

test:
	PYTHONPATH=. $(PYTEST) tests/ -v

smoke:
	PYTHONPATH=. $(PYTEST) tests/test_conference_abstract.py tests/test_pipeline.py -v
	$(PYTHON) scripts/audit_conference_claims.py

reproduce: data validate test audit
	@echo "TCGA-PRAD full scientific pipeline reproduced successfully."

clean:
	rm -rf __pycache__ .pytest_cache
