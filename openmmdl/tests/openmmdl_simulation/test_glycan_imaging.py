"""Tests for the glycan-aware PBC imaging helpers in post_md_conversions."""

import numpy as np
import pytest

# requires MDAnalysis + mdtraj (mdtraj is imported by post_md_conversions)
mda = pytest.importorskip("MDAnalysis")
pytest.importorskip("mdtraj")

from openmmdl.openmmdl_simulation.scripts.post_md_conversions import (  # noqa: E402
    residues_bonded_to_protein,
)


def _build_synthetic_universe():
    """Build a tiny in-memory MDAnalysis Universe consisting of:

    * One Asn residue (atoms N, CA, CB, CG, OD1, ND2)
    * One NAG sugar residue (single C1 atom for compactness)
    * A bond between ND2 and C1 to simulate the N-glycosidic linkage
    * A floating water residue (HOH O) to verify non-bonded non-protein
      atoms are *not* picked up by ``residues_bonded_to_protein``
    """
    n_atoms = 6 + 1 + 1  # Asn + NAG + water
    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=3,
        atom_resindex=[0] * 6 + [1] + [2],
        residue_segindex=[0, 0, 0],
        trajectory=True,
    )
    u.add_TopologyAttr("name", ["N", "CA", "CB", "CG", "OD1", "ND2", "C1", "O"])
    u.add_TopologyAttr("resname", ["ASN", "NAG", "HOH"])
    u.add_TopologyAttr("resid", [1, 2, 3])
    u.add_TopologyAttr("segid", ["A"])

    coords = np.array(
        [
            [0.0, 0.0, 0.0],  # N
            [1.0, 0.0, 0.0],  # CA
            [2.0, 0.0, 0.0],  # CB
            [3.0, 0.0, 0.0],  # CG
            [3.5, 0.7, 0.0],  # OD1
            [3.5, -0.7, 0.0],  # ND2
            [4.7, -1.0, 0.0],  # NAG.C1
            [9.0, 9.0, 9.0],  # HOH.O (far away)
        ],
        dtype=np.float32,
    )
    u.atoms.positions = coords
    # Bond ND2 (index 5) <-> NAG.C1 (index 6)
    u.add_bonds([(5, 6)])
    return u


def test_residues_bonded_to_protein_finds_glycan():
    u = _build_synthetic_universe()
    glycan = residues_bonded_to_protein(u)
    assert glycan.n_atoms == 1
    assert glycan[0].resname == "NAG"


def test_residues_bonded_to_protein_excludes_unbonded_residues():
    u = _build_synthetic_universe()
    glycan = residues_bonded_to_protein(u)
    resnames = {a.resname for a in glycan}
    assert "HOH" not in resnames
    assert "ASN" not in resnames


def test_residues_bonded_to_protein_empty_when_no_bonds():
    # Same atom layout but no bonds at all.
    u = mda.Universe.empty(
        n_atoms=2,
        n_residues=2,
        atom_resindex=[0, 1],
        residue_segindex=[0, 0],
        trajectory=True,
    )
    u.add_TopologyAttr("name", ["ND2", "C1"])
    u.add_TopologyAttr("resname", ["ASN", "NAG"])
    u.add_TopologyAttr("resid", [1, 2])
    u.add_TopologyAttr("segid", ["A"])
    u.atoms.positions = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)
    # No add_bonds call.
    glycan = residues_bonded_to_protein(u)
    assert glycan.n_atoms == 0


def test_residues_bonded_to_protein_handles_no_protein():
    # Universe with only sugar residues -> protein selection is empty.
    u = mda.Universe.empty(
        n_atoms=2,
        n_residues=2,
        atom_resindex=[0, 1],
        residue_segindex=[0, 0],
        trajectory=True,
    )
    u.add_TopologyAttr("name", ["C1", "C1"])
    u.add_TopologyAttr("resname", ["NAG", "NAG"])
    u.add_TopologyAttr("resid", [1, 2])
    u.add_TopologyAttr("segid", ["A"])
    u.atoms.positions = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)
    u.add_bonds([(0, 1)])
    glycan = residues_bonded_to_protein(u)
    assert glycan.n_atoms == 0
