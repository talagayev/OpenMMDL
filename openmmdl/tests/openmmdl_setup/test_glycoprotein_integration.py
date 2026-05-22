"""End-to-end test: glycoprotein PDB upload -> createAmberBashScript output.

Asserts the generated bash script contains the right tleap-mode pieces
(both protein and glycan force fields sourced, rename + bond runtime
helpers emitted, awk filter skipped, etc.). Does not actually execute
tleap.
"""

import tempfile
import textwrap

import pytest

# skip the module if the heavy deps (rdkit/openmm) aren't installed
pytest.importorskip("openmm")
pytest.importorskip("rdkit")
pytest.importorskip("flask")

from openmmdl.openmmdl_setup import openmmdlsetup as M  # noqa: E402

_TINY_GLYCOPROTEIN_PDB = textwrap.dedent("""\
    ATOM      1  N   ASN A 250      10.000  10.000  10.000  1.00  0.00           N
    ATOM      2  CA  ASN A 250      11.000  10.000  10.000  1.00  0.00           C
    ATOM      3  CB  ASN A 250      12.000  10.000  10.000  1.00  0.00           C
    ATOM      4  CG  ASN A 250      13.000  10.000  10.000  1.00  0.00           C
    ATOM      5  OD1 ASN A 250      13.500  10.866  10.000  1.00  0.00           O
    ATOM      6  ND2 ASN A 250      13.500   9.134  10.000  1.00  0.00           N
    HETATM    7  C1  NAG A 600      14.800   8.700  10.000  1.00  0.00           C
    HETATM    8  O4  NAG A 600      15.800   9.700  10.000  1.00  0.00           O
    HETATM    9  N2  NAG A 600      15.500   6.800  10.000  1.00  0.00           N
    CONECT    6    7
    END
    """)


class _FakeSession(dict):
    def get(self, k, d=None):
        return super().get(k, d)


@pytest.fixture
def glycoprot_session(monkeypatch):
    """Set up uploadedFiles + session as if the user had POSTed the form."""
    temp = tempfile.TemporaryFile()
    temp.write(_TINY_GLYCOPROTEIN_PDB.encode())
    temp.seek(0)
    monkeypatch.setattr(M, "uploadedFiles", {"glycoprotFile": [(temp, "tinyglyco.pdb")]})
    fake_session = _FakeSession(
        {
            "rcpType": "glycoprotRcp",
            "glycoprot_prot_ff": "leaprc.protein.ff14SB",
            "glycoprot_glycan_ff": "leaprc.GLYCAM_06j-1",
            "nmLig": False,
            "spLig": False,
            "addType": "addWater",
            "boxType": "cube",
            "dist": "10",
            "water_ff": "tip3p",
            "pos_ion": "Na+",
            "neg_ion": "Cl-",
            "ionConc": "0.15",
        }
    )
    monkeypatch.setattr(M, "session", fake_session)
    yield fake_session


def test_glycoprotein_script_sources_both_force_fields(glycoprot_session):
    script = M.createAmberBashScript()
    assert "source ${rcp_ff}" in script
    assert "source ${glycan_ff}" in script
    assert "prot_ff=leaprc.protein.ff14SB" in script
    assert "glycan_ff=leaprc.GLYCAM_06j-1" in script


def test_glycoprotein_script_uses_awk_cap_filter(glycoprot_session):
    # the awk filter strips ACE/NME cap atoms that tleap rebuilds from templates
    script = M.createAmberBashScript()
    assert "awk" in script
    assert "CH3|HH31|HH32|HH33" in script


def test_glycoprotein_script_uses_pdb4amber_no_reduce(glycoprot_session):
    script = M.createAmberBashScript()
    assert "pdb4amber -i ${rcp_nm}_renamed.pdb" in script
    assert "--no-reduce" in script


def test_glycoprotein_script_emits_rename_and_bond_helpers(glycoprot_session):
    script = M.createAmberBashScript()
    assert "_glyco_rename.py" in script
    assert "_glyco_bonds.py" in script
    assert "tleap.bonds.txt" in script
    assert "write_renamed_pdb" in script
    assert "emit_bond_statements" in script


def test_glycoprotein_script_inlines_renamed_residues(glycoprot_session):
    script = M.createAmberBashScript()
    assert "'0YB'" in script
    assert "rename_map = {" in script


def test_glycoprotein_script_inlines_protein_links(glycoprot_session):
    script = M.createAmberBashScript()
    assert "'ND2'" in script
    assert "'C1'" in script


def test_glycoprotein_script_inject_bonds_after_loadpdb(glycoprot_session):
    script = M.createAmberBashScript()
    # Without a ligand, bonds are inlined into the main tleap.in heredoc.
    assert "$(cat tleap.bonds.txt)" in script
    loadpdb_idx = script.find("system = loadpdb ${rcp_nm}_cnt_rmv.pdb")
    bond_idx = script.find("$(cat tleap.bonds.txt)")
    assert loadpdb_idx != -1 and bond_idx != -1
    assert bond_idx > loadpdb_idx, "bonds must come after loadpdb"


def test_protein_rcp_script_is_unchanged_no_glycoprot_leakage(monkeypatch):
    temp = tempfile.TemporaryFile()
    temp.write(b"ATOM      1  N   ASN A   1      10 10 10  1.00 0.00\nEND\n")
    temp.seek(0)
    monkeypatch.setattr(M, "uploadedFiles", {"protFile": [(temp, "myprot.pdb")]})
    monkeypatch.setattr(
        M,
        "session",
        _FakeSession(
            {
                "rcpType": "protRcp",
                "prot_ff": "leaprc.protein.ff19SB",
                "nmLig": False,
                "spLig": False,
                "addType": "addWater",
                "boxType": "cube",
                "dist": "10",
                "water_ff": "tip3p",
                "pos_ion": "Na+",
                "neg_ion": "Cl-",
                "ionConc": "0.15",
            }
        ),
    )
    script = M.createAmberBashScript()
    # The awk filter must still appear for the existing protRcp branch.
    assert "awk" in script
    # None of the glycoprotein-specific machinery should appear.
    assert "_glyco_rename" not in script
    assert "tleap.bonds.txt" not in script
    assert "glycan_ff" not in script


def test_setAmberOptions_glycoprotRcp_does_not_500_via_flask_client():
    """Regression: previously stored rename_map (with tuple keys) in Flask
    session, which raised TypeError on JSON serialisation -> 500 to the
    user. Drive the full route via the test client and assert HTTP 200.
    """
    import io

    pdb_body = textwrap.dedent("""\
        ATOM      1  N   ASN A 250      10.000  10.000  10.000  1.00  0.00           N
        ATOM      6  ND2 ASN A 250      13.500   9.134  10.000  1.00  0.00           N
        HETATM    7  C1  NAG A 600      14.800   8.700  10.000  1.00  0.00           C
        HETATM    8  O4  NAG A 600      15.800   9.700  10.000  1.00  0.00           O
        CONECT    6    7
        END
        """)

    client = M.app.test_client()
    resp = client.post(
        "/setAmberOptions",
        data={
            "rcpType": "glycoprotRcp",
            "glycoprot_prot_ff": "leaprc.protein.ff14SB",
            "glycoprot_glycan_ff": "leaprc.GLYCAM_06j-1",
            "addType": "addWater",
            "boxType": "cube",
            "dist": "10",
            "water_ff": "tip3p",
            "pos_ion": "Na+",
            "neg_ion": "Cl-",
            "ionConc": "0.15",
            "glycoprotFile": (io.BytesIO(pdb_body.encode()), "tinyglyco.pdb"),
        },
        content_type="multipart/form-data",
    )

    assert resp.status_code == 200, f"Got HTTP {resp.status_code}: {resp.data[:200]!r}"
    body = resp.data.decode("utf-8", "replace")
    assert "Internal Server Error" not in body
    assert "0YB" in body  # detection actually ran
    assert "glycan_ff=leaprc.GLYCAM_06j-1" in body
