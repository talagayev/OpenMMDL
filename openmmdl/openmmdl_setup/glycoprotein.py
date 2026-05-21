"""Glycan detection and GLYCAM renaming for the OpenMMDL AMBER path.

Supports the common N-glycan core sugars (NAG, NDG, BMA, MAN, FUC, FUL).
Residue, linkage and atom-name codes were verified against the
GLYCAM-06j-1 prep/lib files in GLYCAM-Web/gmml; unsupported sugars or
linkage patterns raise NotImplementedError rather than guessing.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass

logger = logging.getLogger(__name__)


# PDB CCD code -> (sugar letter, anomer letter); "A" = alpha, "B" = beta.
SUGAR_TABLE: dict[str, tuple[str, str]] = {
    "NAG": ("Y", "B"),
    "NDG": ("Y", "A"),
    "BMA": ("M", "B"),
    "MAN": ("M", "A"),
    "FUC": ("f", "A"),
    "FUL": ("f", "B"),
}

PDB_SUGAR_RESNAMES: frozenset[str] = frozenset(SUGAR_TABLE.keys())

# Substituted position(s) on a sugar -> GLYCAM linkage code.
LINKAGE_CODE: dict[tuple[int, ...], str] = {
    (): "0",
    (2,): "2",
    (3,): "3",
    (4,): "4",
    (6,): "6",
    (2, 3): "Z",
    (2, 4): "Y",
    (2, 6): "X",
    (3, 4): "W",
    (3, 6): "V",
    (4, 6): "U",
    (2, 3, 4): "T",
    (2, 3, 6): "S",
    (2, 4, 6): "R",
    (3, 4, 6): "Q",
    (2, 3, 4, 6): "P",
}

_HYDROXYL_H_RENAMES = {
    "HO2": "H2O",
    "HO3": "H3O",
    "HO4": "H4O",
    "HO6": "H6O",
}

_NAG_ACETYL_RENAMES = {
    "C7": "C2N",
    "O7": "O2N",
    "C8": "CME",
    "H81": "H1M",
    "H82": "H2M",
    "H83": "H3M",
    "HN2": "H2N",
}

# PDB CCD -> GLYCAM atom-name renames, keyed by original residue name.
SUGAR_ATOM_RENAME: dict[str, dict[str, str]] = {
    "NAG": {**_HYDROXYL_H_RENAMES, **_NAG_ACETYL_RENAMES},
    "NDG": {**_HYDROXYL_H_RENAMES, **_NAG_ACETYL_RENAMES},
    "BMA": dict(_HYDROXYL_H_RENAMES),
    "MAN": dict(_HYDROXYL_H_RENAMES),
    "FUC": dict(_HYDROXYL_H_RENAMES),
    "FUL": dict(_HYDROXYL_H_RENAMES),
}

# A glycosylated Asn is renamed to the GLYCAM NLN template.
ACCEPTOR_RESIDUE_RENAME: dict[str, str] = {
    "ASN": "NLN",
}

ACCEPTOR_ATOM_RENAME: dict[str, dict[str, str]] = {
    "ASN": {"HD2": "HD21"},
}

# HD22 on Asn is replaced by the glycosidic bond and is absent from NLN.
ACCEPTOR_ATOM_DROP: dict[str, frozenset[str]] = {
    "ASN": frozenset({"HD22"}),
}

_PROTEIN_ACCEPTORS = {
    ("ASN", "ND2"),
    ("SER", "OG"),
    ("SEP", "OG"),
    ("THR", "OG1"),
    ("TPO", "OG1"),
    ("TYR", "OH"),
}

PROTEIN_RESNAMES: frozenset[str] = frozenset(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "HID",
        "HIE",
        "HIP",
        "HSD",
        "HSE",
        "HSP",
        "CYX",
        "CYM",
        "ASH",
        "GLH",
        "LYN",
        "SEP",
        "TPO",
    ]
)


@dataclass(frozen=True)
class _Atom:
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resnum: int
    icode: str
    x: float
    y: float
    z: float

    @property
    def res_key(self) -> tuple[str, int, str]:
        return (self.chain, self.resnum, self.icode)


def _parse_pdb(pdb_path: str):
    atoms: list[_Atom] = []
    by_serial: dict[int, _Atom] = {}
    conect: set[tuple[int, int]] = set()

    with open(pdb_path, "r") as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec in ("ATOM", "HETATM"):
                try:
                    serial = int(line[6:11])
                    name = line[12:16].strip()
                    altloc = line[16:17]
                    resname = line[17:20].strip()
                    chain = line[21:22] if len(line) > 21 else " "
                    if not chain.strip():
                        chain = " "
                    resnum = int(line[22:26])
                    icode = line[26:27] if len(line) > 26 else " "
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                a = _Atom(serial, name, altloc, resname, chain, resnum, icode, x, y, z)
                atoms.append(a)
                by_serial[serial] = a
            elif rec == "CONECT":
                try:
                    s = int(line[6:11])
                except ValueError:
                    continue
                for start in (11, 16, 21, 26):
                    chunk = line[start : start + 5].strip() if len(line) > start else ""
                    if not chunk:
                        continue
                    try:
                        t = int(chunk)
                    except ValueError:
                        continue
                    if s != t:
                        conect.add((min(s, t), max(s, t)))

    return atoms, by_serial, sorted(conect)


def _distance(a: _Atom, b: _Atom) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def detect_glycans(pdb_path: str) -> list[dict]:
    """Return the sugar residues found in the PDB.

    Each entry has keys ``chain``, ``resnum``, ``icode``, ``resname`` and
    ``atoms`` (a list of ``(atom_name, serial)`` tuples in file order).
    """
    atoms, _by_serial, _conect = _parse_pdb(pdb_path)

    grouped: dict[tuple[str, int, str], dict] = {}
    order: list[tuple[str, int, str]] = []
    for a in atoms:
        if a.resname not in PDB_SUGAR_RESNAMES:
            continue
        key = a.res_key
        if key not in grouped:
            order.append(key)
            grouped[key] = {
                "chain": a.chain,
                "resnum": a.resnum,
                "icode": a.icode,
                "resname": a.resname,
                "atoms": [],
            }
        grouped[key]["atoms"].append((a.name, a.serial))

    return [grouped[k] for k in order]


def detect_glycosylation_sites(pdb_path: str) -> list[dict]:
    """Find protein-glycan covalent links: CONECT records first, then a
    distance fallback (1.8 A) on acceptor / sugar-C1 atoms."""
    atoms, by_serial, conect = _parse_pdb(pdb_path)

    sites: list[dict] = []
    seen: set = set()

    for s, t in conect:
        if s not in by_serial or t not in by_serial:
            continue
        a, b = by_serial[s], by_serial[t]
        if a.res_key == b.res_key:
            continue
        site = _classify_protein_sugar_bond(a, b)
        if site is None:
            continue
        sig = _site_signature(site)
        if sig not in seen:
            seen.add(sig)
            sites.append(site)

    if sites:
        return sites

    acceptor_atoms = [a for a in atoms if (a.resname, a.name) in _PROTEIN_ACCEPTORS]
    sugar_c1_atoms = [a for a in atoms if a.resname in PDB_SUGAR_RESNAMES and a.name == "C1"]
    cutoff = 1.8
    for acc in acceptor_atoms:
        for c1 in sugar_c1_atoms:
            if acc.res_key == c1.res_key:
                continue
            if _distance(acc, c1) <= cutoff:
                site = _make_site(acc, c1)
                sig = _site_signature(site)
                if sig not in seen:
                    seen.add(sig)
                    sites.append(site)

    return sites


def detect_glycan_glycan_links(pdb_path: str, glycans: list[dict] | None = None) -> list[dict]:
    """Find sugar-sugar covalent links: CONECT records first, then a
    distance fallback (1.7 A) on C1-O atom pairs across sugar residues."""
    del glycans  # accepted for API symmetry
    atoms, by_serial, conect = _parse_pdb(pdb_path)

    links: list[dict] = []
    seen: set = set()

    for s, t in conect:
        if s not in by_serial or t not in by_serial:
            continue
        a, b = by_serial[s], by_serial[t]
        if a.res_key == b.res_key:
            continue
        if a.resname not in PDB_SUGAR_RESNAMES:
            continue
        if b.resname not in PDB_SUGAR_RESNAMES:
            continue
        link = _make_glycan_link(a, b)
        sig = _link_signature(link)
        if sig not in seen:
            seen.add(sig)
            links.append(link)

    if links:
        return links

    c1_atoms = [a for a in atoms if a.resname in PDB_SUGAR_RESNAMES and a.name == "C1"]
    o_atoms = [
        a
        for a in atoms
        if a.resname in PDB_SUGAR_RESNAMES
        and a.name.startswith("O")
        and a.name not in ("O5",)
        and _atom_position(a.name) is not None
    ]
    cutoff = 1.7
    for c1 in c1_atoms:
        for o in o_atoms:
            if c1.res_key == o.res_key:
                continue
            if _distance(c1, o) <= cutoff:
                link = _make_glycan_link(o, c1)
                sig = _link_signature(link)
                if sig not in seen:
                    seen.add(sig)
                    links.append(link)

    return links


def determine_glycam_names(
    glycans: list[dict],
    protein_links: list[dict],
    glycan_glycan_links: list[dict],
) -> dict[tuple[str, int], str]:
    """Map each detected sugar (and glycosylated Asn) to its GLYCAM code.

    Raises NotImplementedError for sugars outside SUGAR_TABLE, non-canonical
    linkages, or substitution patterns absent from LINKAGE_CODE.
    """
    subs: dict[tuple[str, int], set[int]] = defaultdict(set)

    for link in glycan_glycan_links:
        from_atom = link["from_atom"]
        to_atom = link["to_atom"]
        from_key = (link["from_chain"], link["from_resid"])
        to_key = (link["to_chain"], link["to_resid"])

        from_pos = _atom_position(from_atom)
        to_pos = _atom_position(to_atom)

        if from_atom == "C1" and to_pos is not None:
            subs[to_key].add(to_pos)
        elif to_atom == "C1" and from_pos is not None:
            subs[from_key].add(from_pos)
        else:
            raise NotImplementedError(
                f"Glycan-glycan link {link!r} does not match the canonical "
                "C1-to-O sugar linkage pattern. Only standard N-glycan core "
                "linkages (C1 -> O2/O3/O4/O6) are currently supported."
            )

    for site in protein_links:
        if site["sugar_resname"] not in PDB_SUGAR_RESNAMES:
            raise NotImplementedError(
                f"Glycosylation site references unsupported sugar "
                f"{site['sugar_resname']!r}. Supported: "
                f"{sorted(PDB_SUGAR_RESNAMES)}."
            )

    rename: dict[tuple[str, int], str] = {}

    for site in protein_links:
        prot_resname = site["prot_resname"]
        if prot_resname in ACCEPTOR_RESIDUE_RENAME:
            rename[(site["prot_chain"], site["prot_resid"])] = ACCEPTOR_RESIDUE_RENAME[
                prot_resname
            ]
    for entry in glycans:
        resname = entry["resname"]
        if resname not in SUGAR_TABLE:
            raise NotImplementedError(
                f"PDB sugar {resname!r} is not in the supported set "
                f"({sorted(SUGAR_TABLE.keys())}). Extend SUGAR_TABLE only "
                "with codes verified against GLYCAM-Web/gmml."
            )
        sugar_letter, anomer_letter = SUGAR_TABLE[resname]
        key = (entry["chain"], entry["resnum"])
        positions = tuple(sorted(subs.get(key, set())))
        if positions not in LINKAGE_CODE:
            raise NotImplementedError(
                f"Substitution pattern {positions!r} on "
                f"{resname} {entry['chain']}{entry['resnum']} has no "
                "verified GLYCAM linkage code. Edit LINKAGE_CODE only "
                "with entries cross-checked against GLYCAM-Web/gmml."
            )
        rename[key] = f"{LINKAGE_CODE[positions]}{sugar_letter}{anomer_letter}"

    return rename


def write_renamed_pdb(
    input_pdb: str,
    output_pdb: str,
    rename_map: dict[tuple[str, int], str],
) -> None:
    """Write a PDB with glycan / acceptor residues renamed for GLYCAM.

    A ``TER`` is inserted between consecutive different sugar residues so
    tleap does not auto-bond adjacent terminal sugars; the explicit ``bond``
    statements from :func:`emit_bond_statements` reconnect the real
    linkages. ``rename_map`` is the dict from :func:`determine_glycam_names`.
    """
    prev_glycan_key: tuple[str, int] | None = None
    with open(input_pdb, "r") as src, open(output_pdb, "w") as dst:
        for line in src:
            if not line.endswith("\n"):
                line = line + "\n"
            rec = line[0:6].strip()
            if rec in ("ATOM", "HETATM"):
                try:
                    orig_resname = line[17:20].strip()
                    chain = line[21:22]
                    if not chain.strip():
                        chain = " "
                    resnum = int(line[22:26])
                    atom_name = line[12:16].strip()
                except (IndexError, ValueError):
                    dst.write(line)
                    continue
                key = (chain, resnum)
                is_sugar_residue = key in rename_map and orig_resname in SUGAR_ATOM_RENAME

                if is_sugar_residue and prev_glycan_key is not None and prev_glycan_key != key:
                    dst.write("TER\n")

                if key not in rename_map:
                    prev_glycan_key = None
                    dst.write(line)
                    continue

                if atom_name in ACCEPTOR_ATOM_DROP.get(orig_resname, frozenset()):
                    continue

                new_resname = rename_map[key]
                renamed_atom = SUGAR_ATOM_RENAME.get(orig_resname, {}).get(
                    atom_name,
                    ACCEPTOR_ATOM_RENAME.get(orig_resname, {}).get(atom_name, atom_name),
                )

                # rewrite the name (12:16) and resname (17:20) fields
                name_field = _format_atom_name(renamed_atom)
                resname_field = (new_resname + "   ")[:3]
                new_line = line[:12] + name_field + line[16:17] + resname_field + line[20:]
                dst.write(new_line)
                prev_glycan_key = key if is_sugar_residue else None
            elif rec == "TER":
                prev_glycan_key = None
                dst.write(line)
            else:
                dst.write(line)


def emit_bond_statements(
    pdb_path: str,
    protein_links: list[dict],
    glycan_glycan_links: list[dict],
    system_var: str = "system",
) -> list[str]:
    """Return tleap ``bond`` statements for the detected linkages.

    Residue indices are the 1-based sequential order in ``pdb_path``; the
    caller must run tleap on this PDB or a downstream file that preserves
    residue order (e.g. the output of ``pdb4amber``).
    """
    index = _residue_sequential_index(pdb_path)
    statements: list[str] = []

    for site in protein_links:
        p_idx = index.get((site["prot_chain"], site["prot_resid"]))
        s_idx = index.get((site["sugar_chain"], site["sugar_resid"]))
        if p_idx is None or s_idx is None:
            continue
        statements.append(
            f"bond {system_var}.{p_idx}.{site['prot_atom']} {system_var}.{s_idx}.{site['sugar_atom']}"
        )

    for link in glycan_glycan_links:
        f_idx = index.get((link["from_chain"], link["from_resid"]))
        t_idx = index.get((link["to_chain"], link["to_resid"]))
        if f_idx is None or t_idx is None:
            continue
        statements.append(
            f"bond {system_var}.{f_idx}.{link['from_atom']} {system_var}.{t_idx}.{link['to_atom']}"
        )

    return statements


def _classify_protein_sugar_bond(a: _Atom, b: _Atom) -> dict | None:
    if (a.resname, a.name) in _PROTEIN_ACCEPTORS and b.resname in PDB_SUGAR_RESNAMES:
        prot, sugar = a, b
    elif (b.resname, b.name) in _PROTEIN_ACCEPTORS and a.resname in PDB_SUGAR_RESNAMES:
        prot, sugar = b, a
    else:
        return None
    if sugar.name != "C1":
        return None
    return _make_site(prot, sugar)


def _make_site(prot: _Atom, sugar: _Atom) -> dict:
    return {
        "prot_chain": prot.chain,
        "prot_resid": prot.resnum,
        "prot_icode": prot.icode,
        "prot_resname": prot.resname,
        "prot_atom": prot.name,
        "sugar_chain": sugar.chain,
        "sugar_resid": sugar.resnum,
        "sugar_icode": sugar.icode,
        "sugar_resname": sugar.resname,
        "sugar_atom": sugar.name,
    }


def _site_signature(site: dict) -> tuple:
    return (
        site["prot_chain"],
        site["prot_resid"],
        site["prot_atom"],
        site["sugar_chain"],
        site["sugar_resid"],
        site["sugar_atom"],
    )


def _make_glycan_link(a: _Atom, b: _Atom) -> dict:
    return {
        "from_chain": a.chain,
        "from_resid": a.resnum,
        "from_icode": a.icode,
        "from_resname": a.resname,
        "from_atom": a.name,
        "to_chain": b.chain,
        "to_resid": b.resnum,
        "to_icode": b.icode,
        "to_resname": b.resname,
        "to_atom": b.name,
    }


def _link_signature(link: dict) -> tuple:
    return (
        link["from_chain"],
        link["from_resid"],
        link["from_atom"],
        link["to_chain"],
        link["to_resid"],
        link["to_atom"],
    )


def _atom_position(atom_name: str) -> int | None:
    """Substitution position encoded in a non-anomeric ring oxygen name
    (``O3`` -> 3); ``None`` for the ring oxygen O5 and anomeric O1."""
    if len(atom_name) < 2 or atom_name[0] != "O":
        return None
    rest = atom_name[1:]
    if not rest[0].isdigit():
        return None
    try:
        pos = int(rest[0])
    except ValueError:
        return None
    if pos in (2, 3, 4, 6, 7, 8, 9):
        return pos
    return None


def _format_atom_name(name: str) -> str:
    """Format an atom name into the 4-character PDB atom-name field."""
    if len(name) >= 4:
        return name[:4]
    return (" " + name).ljust(4)


def _residue_sequential_index(pdb_path: str) -> dict[tuple[str, int], int]:
    """Map (chain, resnum) -> 1-based sequential index in PDB file order."""
    index: dict[tuple[str, int], int] = {}
    seen_keys: list[tuple[str, int]] = []
    with open(pdb_path, "r") as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                chain = line[21:22]
                if not chain.strip():
                    chain = " "
                resnum = int(line[22:26])
            except (IndexError, ValueError):
                continue
            key = (chain, resnum)
            if key not in index:
                index[key] = len(seen_keys) + 1
                seen_keys.append(key)
    return index
