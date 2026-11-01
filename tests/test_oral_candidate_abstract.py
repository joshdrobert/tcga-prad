import json
import re
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SPRINT_DIR = REPO_ROOT / "research" / "oral_sprint_20261101"
ABSTRACT_DIR = SPRINT_DIR / "submission_package"


def test_oral_candidate_files_exist():
    assert (ABSTRACT_DIR / "TITLE.txt").is_file()
    assert (ABSTRACT_DIR / "BODY.txt").is_file()
    assert (ABSTRACT_DIR / "PORTAL_READY_ABSTRACT.txt").is_file()
    assert (SPRINT_DIR / "CLAIM_LEDGER.csv").is_file()
    assert (SPRINT_DIR / "artifacts" / "ORAL_SPRINT_METRICS.json").is_file()


def test_oral_candidate_character_and_structure_constraints():
    title = (ABSTRACT_DIR / "TITLE.txt").read_text(encoding="utf-8").strip()
    body = (ABSTRACT_DIR / "BODY.txt").read_text(encoding="utf-8").strip()

    title_words = len(title.split())
    assert title_words <= 25, f"Title too long: {title_words} words"

    # Rule 1: Portal body <= 3000 chars with spaces
    body_with_spaces = len(body)
    assert body_with_spaces <= 3000, f"Portal limit exceeded: {body_with_spaces} > 3000"

    # Rule 2: Public call body <= 2500 chars without spaces
    body_no_spaces = len("".join(body.split()))
    assert body_no_spaces <= 2500, f"Public call limit exceeded: {body_no_spaces} > 2500"

    # Total with title without whitespace <= 2500
    total_no_spaces = len("".join((title + body).split()))
    assert total_no_spaces <= 2500, f"Total text without spaces exceeded: {total_no_spaces} > 2500"

    # Required 4 headings
    for heading in ["Aim:", "Material & Methods:", "Results:", "Conclusion:"]:
        assert heading in body, f"Missing required heading: {heading}"

    # Forbidden terms
    forbidden = ["http://", "https://", "@", "Genomic Health", "Decipher", "Oncotype", "Prolaris"]
    for term in forbidden:
        assert term.lower() not in body.lower(), f"Forbidden term found: {term}"


def test_oral_candidate_metrics_exact_agreement():
    body = (ABSTRACT_DIR / "BODY.txt").read_text(encoding="utf-8")
    metrics_path = SPRINT_DIR / "artifacts" / "ORAL_SPRINT_METRICS.json"
    with open(metrics_path) as f:
        m = json.load(f)

    # Cohort sizes and denominators
    assert "187" in body
    assert "314" in body
    assert "19 true recurrence events" in body
    assert "43.2%" in body
    assert "89.5%" in body
    assert "63.2%" in body
    assert "8.16" in body

    # Frozen common-fold metrics
    assert "0.6540" in body
    assert "0.7457" in body
    assert "+0.0917" in body
    assert "0.7069" in body
    assert "0.6917" in body
    assert "-0.0153" in body
    assert "+0.1070" in body

    # Sensitivity comparator
    assert "+0.1168" in body
    assert "-0.0076" in body


def test_final_oral_submission_package_exists_and_conforms():
    pkg_dir = REPO_ROOT / "conference" / "osaka_2027" / "FINAL_ORAL_SUBMISSION_20261101"
    required_files = [
        "FINAL_TITLE.txt",
        "FINAL_ABSTRACT.txt",
        "EXACT_PORTAL_COPY.txt",
        "SCIENTIFIC_EVIDENCE_AUDIT.md",
        "REVIEWER_COMPETITIVENESS_ASSESSMENT.md",
        "ABSTRACT_COMPLIANCE_REPORT.md",
        "AUTHOR_SUBMISSION_CHECKLIST.md",
        "SUBMISSION_DECISION.md",
    ]
    for rf in required_files:
        p = pkg_dir / rf
        assert p.is_file(), f"Missing required final file: {rf}"

    # Verify abstract constraints
    final_abstract = (pkg_dir / "FINAL_ABSTRACT.txt").read_text(encoding="utf-8").strip()
    final_title = (pkg_dir / "FINAL_TITLE.txt").read_text(encoding="utf-8").strip()

    assert len(final_title.split()) <= 25
    assert len(final_abstract) <= 3000
    assert len("".join(final_abstract.split())) <= 2500
    assert len("".join((final_title + final_abstract).split())) <= 2500

    assert "—" not in final_abstract
    assert "--" not in final_abstract
    assert "underscore" not in final_abstract.lower()
    for h in ["Aim:", "Material & Methods:", "Results:", "Conclusion:"]:
        assert h in final_abstract
