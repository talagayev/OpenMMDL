"""Unit tests for openmmdl.openmmdl_setup.glycoprotein."""

import os
import tempfile
import textwrap

import pytest

from openmmdl.openmmdl_setup.glycoprotein import (
    LINKAGE_CODE,
    PDB_SUGAR_RESNAMES,
    SUGAR_ATOM_RENAME,
    SUGAR_TABLE,
    detect_glycan_glycan_links,
    detect_glycans,
    detect_glycosylation_sites,
    determine_glycam_names,
    emit_bond_statements,
    write_renamed_pdb,
)


def _write_pdb(content):
    fd, path = tempfile.mkstemp(suffix=".pdb")
    with os.fdopen(fd, "w") as fh:
        fh.write(content)
    return path


@pytest.fixture
def empty_pdb():
    path = _write_pdb("END\n")
    yield path
    os.unlink(path)


@pytest.fixture
def asn_nag_pdb():
    """Single N-glycan: Asn 250 -- NAG 600 (single sugar, terminal)."""
    pdb = textwrap.dedent("""\
        ATOM      1  N   ASN A 250      10.000  10.000  10.000  1.00  0.00           N
        ATOM      2  CA  ASN A 250      11.000  10.000  10.000  1.00  0.00           C
        ATOM      3  CB  ASN A 250      12.000  10.000  10.000  1.00  0.00           C
        ATOM      4  CG  ASN A 250      13.000  10.000  10.000  1.00  0.00           C
        ATOM      5  OD1 ASN A 250      13.500  10.866  10.000  1.00  0.00           O
        ATOM      6  ND2 ASN A 250      13.500   9.134  10.000  1.00  0.00           N
        HETATM    7  C1  NAG A 600      14.800   8.700  10.000  1.00  0.00           C
        HETATM    8  C2  NAG A 600      15.300   7.700  10.000  1.00  0.00           C
        HETATM    9  O3  NAG A 600      16.300   7.000  10.000  1.00  0.00           O
        HETATM   10  O4  NAG A 600      15.800   9.700  10.000  1.00  0.00           O
        HETATM   11  O6  NAG A 600      14.300  10.700  10.000  1.00  0.00           O
        HETATM   12  N2  NAG A 600      15.500   6.800  10.000  1.00  0.00           N
        HETATM   13  C7  NAG A 600      16.500   6.000  10.000  1.00  0.00           C
        CONECT    6    7
        END
        """)
    path = _write_pdb(pdb)
    yield path
    os.unlink(path)


@pytest.fixture
def n_glycan_core_pdb():
    """Full N-glycan core: Asn -- NAG_1 -- NAG_2 -- BMA -- (MAN_a, MAN_b),
    with a core fucose alpha-1,6 on NAG_1.
    """
    pdb = textwrap.dedent("""\
        ATOM      1  N   ASN A 250      10.000  10.000  10.000  1.00  0.00           N
        ATOM      2  CA  ASN A 250      11.000  10.000  10.000  1.00  0.00           C
        ATOM      3  CB  ASN A 250      12.000  10.000  10.000  1.00  0.00           C
        ATOM      4  CG  ASN A 250      13.000  10.000  10.000  1.00  0.00           C
        ATOM      5  OD1 ASN A 250      13.500  10.866  10.000  1.00  0.00           O
        ATOM      6  ND2 ASN A 250      13.500   9.134  10.000  1.00  0.00           N
        HETATM    7  C1  NAG A 600      14.800   8.700  10.000  1.00  0.00           C
        HETATM    8  O4  NAG A 600      15.800   9.700  10.000  1.00  0.00           O
        HETATM   11  O6  NAG A 600      14.300  10.700  10.000  1.00  0.00           O
        HETATM   16  C1  NAG A 601      17.250   9.700  10.000  1.00  0.00           C
        HETATM   17  O4  NAG A 601      18.250   9.700  10.000  1.00  0.00           O
        HETATM   18  C1  BMA A 602      19.500   9.700  10.000  1.00  0.00           C
        HETATM   19  O3  BMA A 602      20.500   9.700  10.000  1.00  0.00           O
        HETATM   20  O6  BMA A 602      19.500  11.000  10.000  1.00  0.00           O
        HETATM   21  C1  MAN A 603      21.700   9.700  10.000  1.00  0.00           C
        HETATM   22  C1  MAN A 604      19.500  12.200  10.000  1.00  0.00           C
        HETATM   23  C1  FUC A 605      15.500  11.700  10.000  1.00  0.00           C
        CONECT    6    7
        CONECT    8   16
        CONECT   17   18
        CONECT   19   21
        CONECT   20   22
        CONECT   11   23
        END
        """)
    path = _write_pdb(pdb)
    yield path
    os.unlink(path)


def test_sugar_table_anomer_letters_are_a_or_b():
    for resname, (_letter, anomer) in SUGAR_TABLE.items():
        assert anomer in ("A", "B"), f"{resname} has non-A/B anomer letter {anomer!r}"


def test_pdb_sugar_resnames_match_table_keys():
    assert PDB_SUGAR_RESNAMES == frozenset(SUGAR_TABLE.keys())


def test_linkage_codes_for_n_glycan_core():
    assert LINKAGE_CODE[()] == "0"
    assert LINKAGE_CODE[(4,)] == "4"
    assert LINKAGE_CODE[(3,)] == "3"
    assert LINKAGE_CODE[(6,)] == "6"
    assert LINKAGE_CODE[(3, 6)] == "V"
    assert LINKAGE_CODE[(4, 6)] == "U"


def test_nag_atom_rename_table_remaps_acetyl_group():
    rename = SUGAR_ATOM_RENAME["NAG"]
    assert rename["C7"] == "C2N"
    assert rename["O7"] == "O2N"
    assert rename["C8"] == "CME"
    assert rename["HN2"] == "H2N"


def test_mannose_atom_rename_table_only_remaps_hydroxyl_hs():
    rename = SUGAR_ATOM_RENAME["MAN"]
    assert rename["HO3"] == "H3O"
    assert rename["HO4"] == "H4O"
    assert rename["HO6"] == "H6O"
    assert "C7" not in rename
    assert "C1" not in rename


def test_detect_glycans_empty_pdb(empty_pdb):
    assert detect_glycans(empty_pdb) == []


def test_detect_glycans_single_residue(asn_nag_pdb):
    glycans = detect_glycans(asn_nag_pdb)
    assert len(glycans) == 1
    g = glycans[0]
    assert g["resname"] == "NAG"
    assert g["chain"] == "A"
    assert g["resnum"] == 600
    assert g["atoms"][0] == ("C1", 7)


def test_detect_glycans_n_glycan_core(n_glycan_core_pdb):
    glycans = detect_glycans(n_glycan_core_pdb)
    resnames = [g["resname"] for g in glycans]
    assert resnames == ["NAG", "NAG", "BMA", "MAN", "MAN", "FUC"]


def test_detect_glycosylation_sites_conect_based(asn_nag_pdb):
    sites = detect_glycosylation_sites(asn_nag_pdb)
    assert len(sites) == 1
    s = sites[0]
    assert s["prot_resname"] == "ASN"
    assert s["prot_atom"] == "ND2"
    assert s["prot_resid"] == 250
    assert s["sugar_resname"] == "NAG"
    assert s["sugar_atom"] == "C1"
    assert s["sugar_resid"] == 600


def test_detect_glycosylation_sites_distance_fallback():
    # no CONECT records -> distance fallback
    pdb_no_conect = textwrap.dedent("""\
        ATOM      6  ND2 ASN A 250      13.500   9.134  10.000  1.00  0.00           N
        HETATM    7  C1  NAG A 600      14.800   8.700  10.000  1.00  0.00           C
        END
        """)
    path = _write_pdb(pdb_no_conect)
    try:
        sites = detect_glycosylation_sites(path)
        assert len(sites) == 1
        assert sites[0]["prot_resname"] == "ASN"
        assert sites[0]["sugar_resname"] == "NAG"
    finally:
        os.unlink(path)


def test_detect_glycosylation_sites_returns_empty_when_no_acceptor():
    pdb = textwrap.dedent("""\
        HETATM    1  C1  NAG A 600      14.800   8.700  10.000  1.00  0.00           C
        END
        """)
    path = _write_pdb(pdb)
    try:
        assert detect_glycosylation_sites(path) == []
    finally:
        os.unlink(path)


def test_detect_glycan_glycan_links_branched_core(n_glycan_core_pdb):
    links = detect_glycan_glycan_links(n_glycan_core_pdb)
    assert len(links) == 5
    pairs = {
        (
            link["from_resname"] + str(link["from_resid"]),
            link["from_atom"],
            link["to_resname"] + str(link["to_resid"]),
            link["to_atom"],
        )
        for link in links
    }
    assert ("NAG600", "O4", "NAG601", "C1") in pairs
    assert ("NAG600", "O6", "FUC605", "C1") in pairs
    assert ("BMA602", "O3", "MAN603", "C1") in pairs
    assert ("BMA602", "O6", "MAN604", "C1") in pairs


def test_determine_glycam_names_terminal(asn_nag_pdb):
    glycans = detect_glycans(asn_nag_pdb)
    sites = detect_glycosylation_sites(asn_nag_pdb)
    links = detect_glycan_glycan_links(asn_nag_pdb)
    names = determine_glycam_names(glycans, sites, links)
    assert names == {("A", 600): "0YB", ("A", 250): "NLN"}


def test_determine_glycam_names_n_glycan_core(n_glycan_core_pdb):
    glycans = detect_glycans(n_glycan_core_pdb)
    sites = detect_glycosylation_sites(n_glycan_core_pdb)
    links = detect_glycan_glycan_links(n_glycan_core_pdb)
    names = determine_glycam_names(glycans, sites, links)
    assert names[("A", 600)] == "UYB"  # NAG_1: O4 + O6 substituted
    assert names[("A", 601)] == "4YB"  # NAG_2: O4 substituted
    assert names[("A", 602)] == "VMB"  # BMA: O3 + O6 substituted
    assert names[("A", 603)] == "0MA"
    assert names[("A", 604)] == "0MA"
    assert names[("A", 605)] == "0fA"


def test_determine_glycam_names_unsupported_sugar_raises():
    glycans = [{"chain": "A", "resnum": 1, "icode": " ", "resname": "SIA", "atoms": []}]
    with pytest.raises(NotImplementedError) as exc:
        determine_glycam_names(glycans, [], [])
    assert "SIA" in str(exc.value)


def test_write_renamed_pdb_residue_and_atom_renames(asn_nag_pdb):
    out_fd, out_path = tempfile.mkstemp(suffix=".pdb")
    os.close(out_fd)
    try:
        write_renamed_pdb(asn_nag_pdb, out_path, {("A", 600): "0YB"})
        with open(out_path) as f:
            text = f.read()
        assert "0YB" in text
        assert " NAG " not in text
        assert " C2N " in text
        assert " ASN " in text
    finally:
        os.unlink(out_path)


def test_write_renamed_pdb_preserves_conect(asn_nag_pdb):
    out_fd, out_path = tempfile.mkstemp(suffix=".pdb")
    os.close(out_fd)
    try:
        write_renamed_pdb(asn_nag_pdb, out_path, {("A", 600): "0YB"})
        with open(out_path) as f:
            lines = f.readlines()
        conect_lines = [l for l in lines if l.startswith("CONECT")]
        assert len(conect_lines) == 1
        assert "    6    7" in conect_lines[0]
    finally:
        os.unlink(out_path)


def test_emit_bond_statements_uses_sequential_indices(asn_nag_pdb):
    sites = detect_glycosylation_sites(asn_nag_pdb)
    stmts = emit_bond_statements(asn_nag_pdb, sites, [])
    assert stmts == ["bond system.1.ND2 system.2.C1"]


def test_emit_bond_statements_n_glycan_core_count(n_glycan_core_pdb):
    sites = detect_glycosylation_sites(n_glycan_core_pdb)
    links = detect_glycan_glycan_links(n_glycan_core_pdb)
    stmts = emit_bond_statements(n_glycan_core_pdb, sites, links)
    assert len(stmts) == 6  # 1 protein-glycan + 5 glycan-glycan
