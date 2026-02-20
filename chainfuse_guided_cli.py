#!/usr/bin/env python3
"""Guided interactive CLI for selecting receptor/ligand chains step-by-step."""

from __future__ import annotations

import argparse
import re
import shlex
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, Set

import chainfuse_cli as core


def format_non_protein_summary(chain: core.ChainEntry) -> str:
    if not chain.non_protein:
        return "None detected"
    fragments: List[str] = []
    for category, residues in sorted(chain.non_protein.items()):
        label = core.CATEGORY_LABELS.get(category, category.title())
        names = ", ".join(sorted(residues))
        fragments.append(f"{label}: {names}")
    return "; ".join(fragments)


def print_chain_detail(chain: core.ChainEntry) -> None:
    print(f"\n[{chain.master_index}] {chain.file_name} chain {chain.chain_id}")
    print(f"  Nickname: {chain.display_name}")
    print(f"  Residues: {chain.residue_count}")
    print(f"  Non-protein: {format_non_protein_summary(chain)}")


def show_selection_table(
    chains: Sequence[core.ChainEntry],
    selected: Iterable[int],
    disabled: Set[int],
) -> None:
    selected_set = set(selected)
    header = f"{'Idx':>4} {'Role':<7} {'File':<25} {'Chain':<5} {'Residues':>9}  Name"
    print("\nCurrent chain overview:")
    print(header)
    print("-" * len(header))
    for chain in chains:
        badge = " "
        if chain.master_index in selected_set:
            badge = "SELECT"
        elif chain.master_index in disabled:
            badge = "LOCKED"
        print(
            f"[{chain.master_index:>3}] {badge:<7} {chain.file_name:<25} "
            f"{chain.chain_id:<5} {chain.residue_count:>9}  {chain.display_name}"
        )
    print("Legend: SELECT=currently chosen, LOCKED=unavailable for this role.")


def show_selection_summary(role: str, indexes: Sequence[int], chains: Sequence[core.ChainEntry]) -> None:
    if not indexes:
        print(f"\nNo chains chosen for the {role} yet.")
        return
    print(f"\n{role.title()} selection ({len(indexes)} chain(s)):")
    for idx in indexes:
        chain = chains[idx]
        print(
            f"  - [{idx}] {chain.file_name} chain {chain.chain_id} "
            f"({chain.residue_count} residues, non-protein: {format_non_protein_summary(chain)})"
        )


def interactive_chain_choice(
    role: str,
    chains: Sequence[core.ChainEntry],
    disabled: Set[int],
) -> List[int]:
    print(f"\n=== {role.title()} Selection ===")
    print(
        "Type chain indexes (comma/range syntax supported) to toggle them on/off. "
        "Use 'info <idx>' for details, 'list' to reprint the table, and 'done' when finished."
    )
    chosen: List[int] = []
    total = len(chains)
    while True:
        show_selection_table(chains, chosen, disabled)
        show_selection_summary(role, chosen, chains)
        raw = input(f"Choose {role} chains > ").strip()
        lowered = raw.lower()
        if lowered in {"", "list"}:
            continue
        if lowered in {"done", "finish"}:
            if not chosen:
                print("Select at least one chain before finishing.")
                continue
            return chosen
        if lowered.startswith("info"):
            parts = lowered.split()
            if len(parts) != 2:
                print("Usage: info <index>")
                continue
            try:
                idx = int(parts[1])
            except ValueError:
                print("Chain index must be an integer.")
                continue
            if idx < 0 or idx >= total:
                print(f"Index {idx} is out of range.")
                continue
            print_chain_detail(chains[idx])
            continue
        try:
            indexes = core.parse_index_spec(raw, total)
        except ValueError as exc:
            print(f"Error: {exc}")
            continue
        for idx in indexes:
            if idx in disabled:
                print(f"  - Chain {idx} is reserved by another role and cannot be selected.")
                continue
            if idx in chosen:
                chosen.remove(idx)
                print(f"Removed chain {idx} from the {role}.")
            else:
                chosen.append(idx)
                print(f"Added chain {idx} to the {role}.")
                print_chain_detail(chains[idx])


def prompt_choice(prompt: str, options: Sequence[tuple[str, str]], default: str) -> str:
    key_map = {}
    print(f"\n{prompt}")
    for idx, (value, description) in enumerate(options, start=1):
        suffix = " (default)" if value == default else ""
        print(f"  {idx}) {value}: {description}{suffix}")
        key_map[str(idx)] = value
        key_map[value.lower()] = value
    while True:
        raw = input("Select an option (press Enter for default): ").strip().lower()
        if not raw:
            return default
        choice = key_map.get(raw)
        if choice:
            return choice
        print("Invalid selection. Enter the number or name of one of the listed options.")


def describe_plan(role: str, indexes: Sequence[int], mode: str, chains: Sequence[core.ChainEntry]) -> None:
    label = "merged into chain A" if mode == "merge" else "kept as separate chains"
    print(f"- {role.title()}: {len(indexes)} chain(s), {label}")
    for idx in indexes:
        chain = chains[idx]
        print(f"    [{idx}] {chain.display_name} ({chain.file_name} chain {chain.chain_id})")


def prompt_for_pdb_files() -> List[str]:
    cwd = Path.cwd()
    suggestions = sorted(
        entry for entry in cwd.iterdir() if entry.is_file() and entry.suffix.lower() == ".pdb"
    )
    if suggestions:
        print("Detected PDB files in the current directory:")
        for idx, entry in enumerate(suggestions, start=1):
            print(f"  {idx}) {entry.name}")
        print(
            "Type file paths (comma or space separated), "
            "numbers from the list above, or enter 'all' to use every listed file."
        )
        print("Press Enter without typing anything to use the entire detected list.")
    else:
        print("No PDB files detected in the current directory.")
    print("Enter 'quit' to exit the wizard at any time.\n")
    while True:
        raw = input("PDB files > ").strip()
        lowered = raw.lower()
        if not raw:
            if suggestions:
                print("Using all detected PDB files from the current directory.")
                return [str(path.resolve()) for path in suggestions]
            print("Please enter at least one file path.")
            continue
        if lowered in {"quit", "q", "exit"}:
            return []
        if lowered == "all" and suggestions:
            return [str(path.resolve()) for path in suggestions]
        numeric_tokens = [token for token in re.split(r"[,\s]+", raw) if token]
        if (
            suggestions
            and numeric_tokens
            and all(token.isdigit() for token in numeric_tokens)
        ):
            chosen: List[str] = []
            errors = False
            for token in numeric_tokens:
                idx = int(token) - 1
                if idx < 0 or idx >= len(suggestions):
                    print(f"Index {token} is out of range.")
                    errors = True
                    continue
                chosen.append(str(suggestions[idx].resolve()))
            if errors:
                continue
            return chosen
        processed = shlex.split(raw.replace(",", " "))
        if not processed:
            print("Please enter at least one file path.")
            continue
        resolved: List[str] = []
        invalid: List[str] = []
        for token in processed:
            path = Path(token).expanduser()
            if not path.exists() or path.is_dir():
                invalid.append(token)
                continue
            resolved.append(str(path.resolve()))
        if invalid:
            print("Unable to use:", ", ".join(invalid))
            continue
        if resolved:
            return resolved
        print("No valid files detected in that input. Try again.")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Guided wizard for picking receptor/ligand chains interactively.",
    )
    parser.add_argument("pdb_files", nargs="*", help="Input PDB files to inspect")
    parser.add_argument(
        "--output-dir",
        default="chainfuse-output",
        help="Directory for generated receptor/ligand/complex PDB files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting existing files in the output directory",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    print("Welcome to the guided ChainFuse wizard.\n")
    pdb_files = args.pdb_files or prompt_for_pdb_files()
    if not pdb_files:
        print("No PDB files selected. Exiting.")
        return 1

    chains = core.enumerate_chains(pdb_files)
    if not chains:
        print("No chains found in the provided files.")
        return 1

    core.print_chain_summary(chains)
    has_non_protein = core.non_protein_report(chains)

    if has_non_protein:
        print(
            "Some chains include nucleic acids, ligands, ions, or waters. "
            "You can strip them now so the generated structures contain only amino acids."
        )
        strip_non_protein = core.prompt_yes_no("Strip non-protein residues from all outputs?", default=False)
    else:
        strip_non_protein = False
        print("Only protein residues detected; stripping is disabled.\n")

    receptor_indexes = interactive_chain_choice("receptor", chains, disabled=set())
    ligand_indexes = interactive_chain_choice("ligand", chains, disabled=set(receptor_indexes))

    receptor_mode = prompt_choice(
        "How should receptor chains be exported?",
        [
            ("preserve", "Each chain keeps its own identifier"),
            ("merge", "All receptor atoms are merged into a single chain A"),
        ],
        default="preserve",
    )
    ligand_mode = prompt_choice(
        "How should ligand chains be exported?",
        [
            ("preserve", "Each chain keeps its own identifier"),
            ("merge", "All ligand atoms are merged into a single chain A"),
        ],
        default="preserve",
    )
    combined_mode = prompt_choice(
        "How should the receptor+ligand complex be labeled?",
        [
            ("paired", "Receptor becomes chain A and ligand becomes chain B"),
            ("renamed", "Every chain gets renamed sequentially"),
        ],
        default="paired",
    )

    print("\nSummary of your plan:")
    describe_plan("receptor", receptor_indexes, receptor_mode, chains)
    describe_plan("ligand", ligand_indexes, ligand_mode, chains)
    combo_desc = "paired as A/B" if combined_mode == "paired" else "renamed sequentially"
    print(f"- Combined complex will be {combo_desc}.")
    if strip_non_protein:
        print("- Non-protein residues will be stripped from every output.")
    else:
        print("- Non-protein residues will be kept as they appeared in the source files.")

    if not core.prompt_yes_no("\nProceed with generating the PDB files?", default=True):
        print("Aborted before writing files.")
        return 0

    output_dir = Path(args.output_dir)
    receptor_path = output_dir / "receptor.pdb"
    ligand_path = output_dir / "ligand.pdb"
    combined_path = output_dir / "complex.pdb"

    print("\nBuilding receptor structure...")
    receptor_atoms = core.build_component_structure(
        [chains[i] for i in receptor_indexes],
        merge=receptor_mode == "merge",
        merged_label="A",
        strip_non_protein=strip_non_protein,
    )
    if not receptor_atoms:
        raise ValueError("Receptor selection produced no atoms. Try different chains.")

    print("Building ligand structure...")
    ligand_atoms = core.build_component_structure(
        [chains[i] for i in ligand_indexes],
        merge=ligand_mode == "merge",
        merged_label="A",
        strip_non_protein=strip_non_protein,
    )
    if not ligand_atoms:
        raise ValueError("Ligand selection produced no atoms. Try different chains.")

    print("Building combined complex...")
    combined_atoms = core.build_combined_structure(
        [chains[i] for i in receptor_indexes],
        [chains[i] for i in ligand_indexes],
        combined_mode,
        strip_non_protein,
    )

    core.ensure_output(receptor_path, args.overwrite)
    core.ensure_output(ligand_path, args.overwrite)
    core.ensure_output(combined_path, args.overwrite)

    core.write_pdb(receptor_path, receptor_atoms)
    core.write_pdb(ligand_path, ligand_atoms)
    core.write_pdb(combined_path, combined_atoms)

    print("\nDone! Generated files:")
    print(f"  Receptor: {receptor_path}")
    print(f"  Ligand:   {ligand_path}")
    print(f"  Complex:  {combined_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
