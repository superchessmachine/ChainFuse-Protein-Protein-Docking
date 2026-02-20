#!/usr/bin/env python3
"""Interactive CLI tool for inspecting and remixing PDB chains."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple


AMINO_ACIDS = {
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
    "SEC",
    "PYL",
}


@dataclass
class AtomRecord:
    record_name: str
    serial: int
    atom_name: str
    alt_loc: str
    res_name: str
    chain_id: str
    res_seq: int
    ins_code: str
    x: float
    y: float
    z: float
    occupancy: float
    temp_factor: float
    element: str
    charge: str


CATEGORY_LABELS = {
    "nucleic": "DNA/RNA",
    "ion": "Ions",
    "water": "Water",
    "ligand": "Ligands",
}


NUCLEIC_RESIDUES = {
    "A",
    "DA",
    "ADE",
    "C",
    "DC",
    "CYT",
    "G",
    "DG",
    "GUA",
    "T",
    "DT",
    "THY",
    "U",
    "DU",
    "URA",
    "I",
    "DI",
    "PSU",
}


ION_RESIDUES = {
    "MG",
    "MN",
    "ZN",
    "CA",
    "NA",
    "K",
    "CL",
    "FE",
    "CU",
    "CO",
    "NI",
    "BR",
    "I",
    "IOD",
}


WATER_RESIDUES = {"HOH", "H2O", "WAT"}


@dataclass
class ChainEntry:
    master_index: int
    file_path: Path
    file_name: str
    chain_id: str
    display_name: str
    residue_count: int
    atoms: List[AtomRecord]
    non_protein: Dict[str, set[str]]


def parse_atom_line(line: str) -> Optional[AtomRecord]:
    record = line[0:6].strip()
    if record not in {"ATOM", "HETATM"}:
        return None
    try:
        serial = int(line[6:11])
    except ValueError:
        serial = 0
    atom_name = line[12:16].strip()
    alt_loc = line[16:17].strip()
    res_name = line[17:20].strip()
    chain_id = line[21:22].strip() or "?"
    try:
        res_seq = int(line[22:26])
    except ValueError:
        res_seq = 0
    ins_code = line[26:27].strip()
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        raise ValueError(f"Malformed coordinates in line: {line.rstrip()}\n")
    try:
        occupancy = float(line[54:60])
    except ValueError:
        occupancy = 1.0
    try:
        temp_factor = float(line[60:66])
    except ValueError:
        temp_factor = 0.0
    element = line[76:78].strip()
    charge = line[78:80].strip()
    return AtomRecord(
        record_name=record,
        serial=serial,
        atom_name=atom_name,
        alt_loc=alt_loc,
        res_name=res_name,
        chain_id=chain_id,
        res_seq=res_seq,
        ins_code=ins_code,
        x=x,
        y=y,
        z=z,
        occupancy=occupancy,
        temp_factor=temp_factor,
        element=element,
        charge=charge,
    )


def classify_non_protein(res_name: str) -> str:
    name = res_name.strip().upper()
    if not name:
        return "ligand"
    if name in WATER_RESIDUES:
        return "water"
    if name in NUCLEIC_RESIDUES or (len(name) == 2 and name[0] in {"D", "R"} and name[1] in "ACGUIT"):
        return "nucleic"
    if name in ION_RESIDUES or (len(name) <= 2 and name.isalpha()):
        return "ion"
    return "ligand"


def parse_compnd_names(lines: Sequence[str]) -> Dict[str, str]:
    compnd_lines = [line[10:80].rstrip() for line in lines if line.startswith("COMPND")]
    if not compnd_lines:
        return {}
    text = " ".join(compnd_lines)
    blocks = re.split(r"MOL_ID:\s*", text, flags=re.IGNORECASE)
    chain_names: Dict[str, str] = {}
    for block in blocks[1:]:
        block = block.strip()
        if not block:
            continue
        match = re.match(r"(\d+)\s*;?(.*)", block, re.IGNORECASE | re.DOTALL)
        if not match:
            continue
        mol_id, tail = match.groups()
        entries = [segment.strip() for segment in tail.split(";") if segment.strip()]
        fields: Dict[str, str] = {}
        for entry in entries:
            if ":" in entry:
                key, value = entry.split(":", 1)
                fields[key.strip().upper()] = value.strip()
        name = fields.get("MOLECULE") or fields.get("SYNONYM") or f"Molecule {mol_id}"
        chains_field = fields.get("CHAIN")
        if not chains_field:
            continue
        for chain in chains_field.split(","):
            chain_id = chain.strip()
            if chain_id:
                chain_names[chain_id] = name
    return chain_names


def load_pdb(path: Path) -> List[ChainEntry]:
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    chain_names = parse_compnd_names(lines)
    atoms_by_chain: Dict[str, List[AtomRecord]] = {}
    aa_residues: Dict[str, set[Tuple[str, int, str]]] = {}
    residues_all: Dict[str, set[Tuple[str, int, str]]] = {}
    non_protein: Dict[str, Dict[str, set[str]]] = {}
    for line in lines:
        atom = parse_atom_line(line)
        if not atom:
            continue
        atoms_by_chain.setdefault(atom.chain_id, []).append(atom)
        key = (atom.res_name, atom.res_seq, atom.ins_code)
        residues_all.setdefault(atom.chain_id, set()).add(key)
        if atom.res_name.upper() in AMINO_ACIDS:
            aa_residues.setdefault(atom.chain_id, set()).add(key)
        else:
            category = classify_non_protein(atom.res_name)
            bucket = non_protein.setdefault(atom.chain_id, {})
            bucket.setdefault(category, set()).add(atom.res_name.upper())
    entries: List[ChainEntry] = []
    for chain_id, atoms in atoms_by_chain.items():
        display_name = chain_names.get(chain_id) or f"{path.name} chain {chain_id}"
        residue_count = len(aa_residues.get(chain_id, set())) or len(residues_all.get(chain_id, set()))
        entries.append(
            ChainEntry(
                master_index=-1,
                file_path=path,
                file_name=path.name,
                chain_id=chain_id,
                display_name=display_name,
                residue_count=residue_count,
                atoms=atoms,
                non_protein=non_protein.get(chain_id, {}),
            )
        )
    return entries


def enumerate_chains(pdb_files: Sequence[str]) -> List[ChainEntry]:
    chains: List[ChainEntry] = []
    for pdb_file in pdb_files:
        path = Path(pdb_file)
        if not path.exists():
            raise FileNotFoundError(f"Missing PDB file: {pdb_file}")
        for chain in load_pdb(path):
            chain.master_index = len(chains)
            chains.append(chain)
    return chains


def print_chain_summary(chains: Sequence[ChainEntry]) -> None:
    print("Loaded chains:\n")
    header = f"{'Idx':>4}  {'File':<25} {'Chain':<5} {'Residues':>9}  Name"
    print(header)
    print("-" * len(header))
    for chain in chains:
        print(
            f"[{chain.master_index:>3}]  {chain.file_name:<25} {chain.chain_id:<5} {chain.residue_count:>9}  {chain.display_name}"
        )
    print()


def non_protein_report(chains: Sequence[ChainEntry]) -> bool:
    has_non_protein = False
    lines: List[str] = []
    for chain in chains:
        if not chain.non_protein:
            continue
        has_non_protein = True
        fragments = []
        for category, residues in sorted(chain.non_protein.items()):
            names = ", ".join(sorted(residues))
            label = CATEGORY_LABELS.get(category, category.title())
            fragments.append(f"{label}: {names}")
        joined = "; ".join(fragments)
        lines.append(
            f"  - [{chain.master_index}] {chain.file_name} chain {chain.chain_id}: {joined}"
        )
    if has_non_protein:
        print("Non-protein residues detected:")
        for entry in lines:
            print(entry)
        print()
    else:
        print("No non-protein residues (DNA/ligands/ions) detected.\n")
    return has_non_protein


def parse_index_spec(spec: str, total: int) -> List[int]:
    items: List[int] = []
    tokens = [token for token in re.split(r"[\s,]+", spec.strip()) if token]
    for token in tokens:
        if "-" in token:
            start_str, end_str = token.split("-", 1)
            start = int(start_str)
            end = int(end_str)
            if start > end:
                start, end = end, start
            items.extend(range(start, end + 1))
        else:
            items.append(int(token))
    cleaned: List[int] = []
    seen = set()
    for idx in items:
        if idx < 0 or idx >= total:
            raise ValueError(f"Index {idx} out of range (0-{total - 1})")
        if idx not in seen:
            cleaned.append(idx)
            seen.add(idx)
    return cleaned


def prompt_for_indexes(prompt: str, total: int) -> List[int]:
    while True:
        try:
            raw = input(prompt).strip()
            if not raw:
                print("Please enter at least one index.")
                continue
            return parse_index_spec(raw, total)
        except ValueError as exc:
            print(f"Error: {exc}")


def prompt_yes_no(prompt: str, default: bool = False) -> bool:
    suffix = " [Y/n] " if default else " [y/N] "
    while True:
        raw = input(prompt + suffix).strip().lower()
        if not raw:
            return default
        if raw in {"y", "yes"}:
            return True
        if raw in {"n", "no"}:
            return False
        print("Please respond with 'y' or 'n'.")


def chain_id_stream() -> Iterator[str]:
    symbols = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    for symbol in symbols:
        yield symbol
    raise RuntimeError("Ran out of unique chain labels to assign.")


def filter_chain_atoms(chain: ChainEntry, strip_non_protein: bool) -> List[AtomRecord]:
    if not strip_non_protein:
        return chain.atoms
    return [atom for atom in chain.atoms if atom.res_name.upper() in AMINO_ACIDS]


def remap_chain_atoms(
    atoms: Sequence[AtomRecord],
    new_chain_id: str,
    start_res_seq: int = 1,
) -> Tuple[List[AtomRecord], int]:
    residue_map: Dict[Tuple[str, int, str], int] = {}
    current = start_res_seq
    remapped: List[AtomRecord] = []
    for atom in atoms:
        key = (atom.res_name, atom.res_seq, atom.ins_code)
        if key not in residue_map:
            residue_map[key] = current
            current += 1
        remapped.append(
            replace(
                atom,
                chain_id=new_chain_id,
                res_seq=residue_map[key],
                ins_code="",
            )
        )
    return remapped, current


def merge_chains(chains: Sequence[Sequence[AtomRecord]], target_chain_id: str) -> List[AtomRecord]:
    aggregated: List[AtomRecord] = []
    next_residue = 1
    for atoms in chains:
        if not atoms:
            continue
        remapped, next_residue = remap_chain_atoms(atoms, target_chain_id, start_res_seq=next_residue)
        aggregated.extend(remapped)
    return aggregated


def assign_unique_chain_atoms(chains: Sequence[Sequence[AtomRecord]], generator: Iterator[str] | None = None) -> List[AtomRecord]:
    generator = generator or chain_id_stream()
    aggregated: List[AtomRecord] = []
    for atoms in chains:
        if not atoms:
            continue
        try:
            new_chain = next(generator)
        except StopIteration:
            raise RuntimeError("Unable to assign additional chain labels.") from None
        remapped, _ = remap_chain_atoms(atoms, new_chain, start_res_seq=1)
        aggregated.extend(remapped)
    return aggregated


def format_atom_line(atom: AtomRecord, serial: int) -> str:
    atom_name = atom.atom_name
    if len(atom_name) < 4:
        atom_name = atom_name.rjust(4)
    else:
        atom_name = atom_name[:4]
    chain_id = (atom.chain_id or "?")[:1]
    ins_code = (atom.ins_code or " ")[:1]
    element = (atom.element or "").rjust(2)
    charge = (atom.charge or "").rjust(2)
    return (
        f"{atom.record_name:<6}{serial:5d} {atom_name}{atom.alt_loc or ' '}{atom.res_name:>3} "
        f"{chain_id}{atom.res_seq:4d}{ins_code}   {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{atom.occupancy:6.2f}{atom.temp_factor:6.2f}          {element}{charge}"
    )


def write_pdb(path: Path, atoms: Sequence[AtomRecord]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for serial, atom in enumerate(atoms, start=1):
            handle.write(format_atom_line(atom, serial) + "\n")
        handle.write("END\n")


def build_component_structure(
    chains: Sequence[ChainEntry],
    *,
    merge: bool,
    merged_label: str,
    strip_non_protein: bool,
) -> List[AtomRecord]:
    if not chains:
        return []
    atom_groups = [filter_chain_atoms(chain, strip_non_protein) for chain in chains]
    if strip_non_protein:
        emptied = [chain for chain, atoms in zip(chains, atom_groups) if not atoms]
        if emptied:
            labels = ", ".join(f"[{c.master_index}] {c.chain_id}" for c in emptied)
            print(
                f"Warning: chains {labels} contained no protein residues and will be skipped."
            )
    if merge:
        return merge_chains(atom_groups, merged_label)
    return assign_unique_chain_atoms(atom_groups)


def build_combined_structure(
    receptor: Sequence[ChainEntry],
    ligand: Sequence[ChainEntry],
    mode: str,
    strip_non_protein: bool,
) -> List[AtomRecord]:
    receptor_atoms = [filter_chain_atoms(chain, strip_non_protein) for chain in receptor]
    ligand_atoms = [filter_chain_atoms(chain, strip_non_protein) for chain in ligand]
    if mode == "paired":
        atoms: List[AtomRecord] = []
        atoms.extend(merge_chains(receptor_atoms, "A"))
        atoms.extend(merge_chains(ligand_atoms, "B"))
        return atoms
    if mode == "renamed":
        generator = chain_id_stream()
        atoms = assign_unique_chain_atoms(receptor_atoms, generator)
        atoms.extend(assign_unique_chain_atoms(ligand_atoms, generator))
        return atoms
    raise ValueError(f"Unsupported combined mode: {mode}")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect PDB chains and build receptor/ligand/complex structures.",
    )
    parser.add_argument("pdb_files", nargs="+", help="Input PDB files")
    parser.add_argument(
        "--receptor",
        help="Comma or dash separated list of master indexes to mark as receptor",
    )
    parser.add_argument(
        "--ligand",
        help="Comma or dash separated list of master indexes to mark as ligand",
    )
    parser.add_argument(
        "--receptor-mode",
        choices=["preserve", "merge"],
        default="preserve",
        help="Preserve separate chains or merge receptor into one chain",
    )
    parser.add_argument(
        "--ligand-mode",
        choices=["preserve", "merge"],
        default="preserve",
        help="Preserve separate chains or merge ligand into one chain",
    )
    parser.add_argument(
        "--combined-mode",
        choices=["paired", "renamed"],
        default="paired",
        help="How to label chains in the combined structure",
    )
    parser.add_argument(
        "--strip-non-protein",
        dest="strip_non_protein",
        action="store_true",
        help="Remove non-protein residues (DNA/ligands/ions) without prompting",
    )
    parser.add_argument(
        "--keep-non-protein",
        dest="strip_non_protein",
        action="store_false",
        help="Always keep non-protein residues without prompting",
    )
    parser.set_defaults(strip_non_protein=None)
    parser.add_argument(
        "--output-dir",
        default="chainfuse-output",
        help="Where to place generated PDB files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting existing output files",
    )
    return parser.parse_args(argv)


def ensure_output(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    chains = enumerate_chains(args.pdb_files)
    if not chains:
        print("No chains found in the provided files.", file=sys.stderr)
        return 1
    print_chain_summary(chains)
    has_non_protein = non_protein_report(chains)
    strip_non_protein = args.strip_non_protein
    if strip_non_protein is None:
        if has_non_protein:
            strip_non_protein = prompt_yes_no(
                "Remove non-protein residues (DNA/ligands/ions)?",
                default=False,
            )
        else:
            strip_non_protein = False
    if strip_non_protein:
        print("Will remove non-protein residues from all generated structures.\n")
    elif has_non_protein:
        print("Keeping non-protein residues in generated structures.\n")
    total = len(chains)
    if args.receptor:
        receptor_indexes = parse_index_spec(args.receptor, total)
    else:
        receptor_indexes = prompt_for_indexes("Enter receptor chain indexes: ", total)
    if args.ligand:
        ligand_indexes = parse_index_spec(args.ligand, total)
    else:
        ligand_indexes = prompt_for_indexes("Enter ligand chain indexes: ", total)
    overlap = set(receptor_indexes) & set(ligand_indexes)
    if overlap:
        raise ValueError(f"Receptor and ligand selections overlap: {sorted(overlap)}")
    receptor_chains = [chains[i] for i in receptor_indexes]
    ligand_chains = [chains[i] for i in ligand_indexes]
    if not receptor_chains:
        raise ValueError("Receptor selection is empty.")
    if not ligand_chains:
        raise ValueError("Ligand selection is empty.")

    output_dir = Path(args.output_dir)
    receptor_path = output_dir / "receptor.pdb"
    ligand_path = output_dir / "ligand.pdb"
    combined_path = output_dir / "complex.pdb"

    receptor_atoms = build_component_structure(
        receptor_chains,
        merge=args.receptor_mode == "merge",
        merged_label="A",
        strip_non_protein=strip_non_protein,
    )
    if not receptor_atoms:
        raise ValueError("Receptor selection contains no protein atoms after filtering.")
    ligand_atoms = build_component_structure(
        ligand_chains,
        merge=args.ligand_mode == "merge",
        merged_label="A",
        strip_non_protein=strip_non_protein,
    )
    if not ligand_atoms:
        raise ValueError("Ligand selection contains no protein atoms after filtering.")
    combined_atoms = build_combined_structure(
        receptor_chains,
        ligand_chains,
        args.combined_mode,
        strip_non_protein,
    )
    if not combined_atoms:
        raise ValueError("Combined structure contains no atoms after filtering selections.")

    ensure_output(receptor_path, args.overwrite)
    ensure_output(ligand_path, args.overwrite)
    ensure_output(combined_path, args.overwrite)

    write_pdb(receptor_path, receptor_atoms)
    write_pdb(ligand_path, ligand_atoms)
    write_pdb(combined_path, combined_atoms)

    print("Generated files:")
    print(f"  Receptor: {receptor_path}")
    print(f"  Ligand:   {ligand_path}")
    print(f"  Complex:  {combined_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
