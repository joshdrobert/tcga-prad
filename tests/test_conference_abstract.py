import re
from pathlib import Path

def test_osaka_2027_conference_package_exists():
    base_dir = Path(__file__).resolve().parent.parent / "conference" / "osaka_2027"
    assert (base_dir / "READY_TO_PASTE_ABSTRACT.txt").is_file()
    assert (base_dir / "SUBMISSION_PACKET.md").is_file()
    assert (base_dir / "SCIENTIFIC_FIT_REVIEW.md").is_file()
    assert (base_dir / "PORTAL_CHECKLIST.md").is_file()

def test_osaka_2027_abstract_character_and_structure_constraints():
    abstract_path = Path(__file__).resolve().parent.parent / "conference" / "osaka_2027" / "READY_TO_PASTE_ABSTRACT.txt"
    content = abstract_path.read_text(encoding="utf-8")

    # Extract abstract body
    assert "ABSTRACT TEXT (PASTE BELOW INTO PORTAL):" in content
    parts = content.split("ABSTRACT TEXT (PASTE BELOW INTO PORTAL):")
    abstract_body = parts[1].strip()

    # Verify mandatory headings
    assert "Aim:" in abstract_body
    assert "Material & Methods:" in abstract_body
    assert "Results:" in abstract_body
    assert "Conclusion:" in abstract_body

    # Check character count excluding whitespace
    chars_no_spaces = len(abstract_body.replace(" ", "").replace("\n", "").replace("\r", "").replace("\t", ""))
    assert chars_no_spaces <= 2500, f"Body character count excluding spaces ({chars_no_spaces}) exceeds limit of 2,500"

    # Check title extraction & constraints
    title_match = re.search(r"TITLE:\s*\n([^\n]+)", content)
    assert title_match is not None
    title = title_match.group(1).strip()
    total_with_title_no_spaces = chars_no_spaces + len(title.replace(" ", ""))
    assert total_with_title_no_spaces <= 2500, f"Full abstract with title ({total_with_title_no_spaces}) exceeds 2,500 characters"

def test_osaka_2027_scientific_integrity_and_data_anchors():
    abstract_path = Path(__file__).resolve().parent.parent / "conference" / "osaka_2027" / "READY_TO_PASTE_ABSTRACT.txt"
    body = abstract_path.read_text(encoding="utf-8")

    # Required data points
    assert "314" in body
    assert "261" in body
    assert "187" in body
    assert "0.7104" in body
    assert "0.7028" in body
    assert "-0.0076" in body
    assert "0.1673" in body
    assert "0.1672" in body
    assert "0.756" in body
    assert "0.511" in body
    assert "248" in body
    assert "511.0" in body
    assert "29" in body

    # Prohibited hallucinated phrases
    assert "Antigravity Biomedical Engineering Consortium" not in body
    assert "TCIA" not in body
    assert "PSMA-PET segmentation" not in body
