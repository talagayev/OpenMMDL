import sys
import types

import pytest

pd = pytest.importorskip("pandas")


# Provide light-weight stand-ins for optional heavy dependencies so the module
# can be imported in environments where MDAnalysis/ProLIF are unavailable.
if "MDAnalysis" not in sys.modules:
    mda_stub = types.ModuleType("MDAnalysis")

    class _Universe:  # pragma: no cover - simple import shim
        def __init__(self, *args, **kwargs):
            raise NotImplementedError

    def _merge(*args, **kwargs):  # pragma: no cover - simple import shim
        raise NotImplementedError

    mda_stub.Universe = _Universe
    mda_stub.Merge = _merge
    mda_stub.analysis = types.SimpleNamespace(
        rms=types.SimpleNamespace(RMSD=None),
    )
    sys.modules["MDAnalysis"] = mda_stub

if "prolif" not in sys.modules:
    prolif_stub = types.ModuleType("prolif")

    class _Fingerprint:  # pragma: no cover - simple import shim
        def __init__(self, *args, **kwargs):
            pass

        def run(self, *args, **kwargs):
            raise NotImplementedError

        def to_dataframe(self, *args, **kwargs):
            raise NotImplementedError

    prolif_stub.Fingerprint = _Fingerprint
    sys.modules["prolif"] = prolif_stub


from openmmdl.openmmdl_analysis.analysis.interactions import (
    InteractionAnalyzer,
    _AtomInfo,
    _create_empty_row,
    _format_atom_indices,
    _parse_residue_label,
)


def test_parse_residue_label_handles_common_formats():
    chain, resid, resname = _parse_residue_label("A:LEU:132")
    assert chain == "A"
    assert resid == 132
    assert resname == "LEU"

    chain, resid, resname = _parse_residue_label("A:LEU132")
    assert chain == "A"
    assert resid is None
    assert "LEU" in resname

    chain, resid, resname = _parse_residue_label("LIG1")
    assert chain is not None
    assert resid is None
    assert resname

    chain, resid, resname = _parse_residue_label(None)
    assert chain is None
    assert resid is None
    assert resname is None


def test_format_atom_indices_returns_unique_and_sorted():
    atoms = [
        _AtomInfo(index=7, resid=None, resname=None, chain=None, name=None),
        _AtomInfo(index=2, resid=None, resname=None, chain=None, name=None),
        _AtomInfo(index=7, resid=None, resname=None, chain=None, name=None),
    ]
    assert _format_atom_indices(atoms) == "2,7"

    single_atom = [_AtomInfo(index=5, resid=None, resname=None, chain=None, name=None)]
    assert _format_atom_indices(single_atom) == 5

    assert _format_atom_indices([]) is None


@pytest.fixture
def analyzer():
    instance = InteractionAnalyzer.__new__(InteractionAnalyzer)
    instance.special = None
    return instance


def test_convert_prolif_dataframe_translates_matches(analyzer):
    frame = 5
    acceptor_entry = {
        "protein_atoms": [
            {"index": 5, "resid": 132, "resname": "TYR", "chain": "A", "name": "OH"}
        ],
        "ligand_atoms": [
            {"index": 101, "resid": 1, "resname": "LIG", "chain": "X", "name": "O1"}
        ],
        "metadata": {},
    }
    donor_entry = {
        "protein_atoms": [
            {"index": 7, "resid": 133, "resname": "HIS", "chain": "A", "name": "NH"}
        ],
        "ligand_atoms": [
            {"index": 102, "resid": 1, "resname": "LIG", "chain": "X", "name": "H1"}
        ],
        "metadata": {},
    }

    columns = pd.MultiIndex.from_tuples(
        [
            ("HBAcceptor", "A:TYR132", "LIG1"),
            ("HBDonor", "A:HIS133", "LIG1"),
        ]
    )
    prolif_df = pd.DataFrame(
        {
            columns[0]: [[acceptor_entry]],
            columns[1]: [[donor_entry]],
        },
        index=[frame],
    )

    result = analyzer._convert_prolif_dataframe(prolif_df, frame)

    assert len(result) == 2
    assert set(result["INTERACTION"]) == {"hbond"}

    rows = result.set_index("Prot_partner")

    acceptor_row = rows.loc["132TYRA"]
    assert acceptor_row["FRAME"] == frame
    assert acceptor_row["PROTISDON"] is False
    assert acceptor_row["DONORIDX"] == 101
    assert acceptor_row["ACCEPTORIDX"] == 5
    assert acceptor_row["RESCHAIN"] == "A"
    assert acceptor_row["RESNR_LIG"] == 1
    assert acceptor_row["RESTYPE_LIG"] == "LIG"
    assert acceptor_row["RESCHAIN_LIG"] == "X"

    donor_row = rows.loc["133HISA"]
    assert donor_row["FRAME"] == frame
    assert donor_row["PROTISDON"] is True
    assert donor_row["DONORIDX"] == 7
    assert donor_row["ACCEPTORIDX"] == 102


def test_convert_prolif_dataframe_returns_empty_template(analyzer):
    result = analyzer._convert_prolif_dataframe(pd.DataFrame(), frame=0)
    assert result.empty
    assert list(result.columns) == list(_create_empty_row().keys())
