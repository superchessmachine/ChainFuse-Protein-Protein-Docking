# ChainFuse-Protein-Protein-Docking

`chainfuse_cli.py` is an interactive CLI for inspecting the chains in one or
more PDB files, selecting which chains belong to the receptor and ligand, and
writing out trimmed PDB files (receptor, ligand, and a combined complex)
according to the desired merging/renaming strategy.

`chainfuse_guided_cli.py` offers the same output capabilities but walks you
through every decision with a CLI wizard. It previews each chain as you
toggle it on/off for the receptor or ligand and prints a running summary so
you always know how your complex will be constructed. Run it with zero
arguments to be prompted for the PDB files interactively.

## Requirements

- Python 3.9+

## Usage

```
python chainfuse_cli.py Bre1-Lge1-helical-bundle.pdb EC-noDNARNA.pdb
```

When executed, the tool prints a table containing every detected chain across
the provided PDB files. Each chain is assigned a master index which you can use
to interactively select receptor and ligand chains. After choosing the chain
indexes, the program produces three PDB files in `chainfuse-output/` by default:

- `receptor.pdb` – all selected receptor chains (kept separate or merged)
- `ligand.pdb` – all selected ligand chains (kept separate or merged)
- `complex.pdb` – receptor and ligand combined, with either chains A/B or
  individual renamed chains

During execution the CLI scans every chain for non-protein residues (DNA, ions,
ligands, etc.), reports what it finds, and asks whether you want to strip those
atoms before generating the receptor/ligand/complex PDBs.

### Guided wizard

```
python chainfuse_guided_cli.py
```

The wizard will first let you choose PDB inputs (auto-detecting any in the
current directory if you just press Enter) and then walk you through every
decision. It prints a live chain table, lets you inspect each chain in detail,
and highlights any residues that are not part of the protein backbone. Chains
can be toggled interactively and the receptor/ligand/complex build plan is
summarized before any files are written.

### Options

```
python chainfuse_cli.py [OPTIONS] <pdb files>

Options:
  --receptor <indexes>     Pre-select receptor chains (comma/dash format)
  --ligand <indexes>       Pre-select ligand chains (comma/dash format)
  --receptor-mode <mode>   'preserve' (default) or 'merge' chains in receptor
  --ligand-mode <mode>     'preserve' or 'merge' chains in ligand
  --combined-mode <mode>   'paired' (receptor=A, ligand=B) or 'renamed'
  --strip-non-protein      Remove DNA/ligand/ion residues without prompting
  --keep-non-protein       Always keep non-protein residues without prompting
  --output-dir <path>      Output directory (default: chainfuse-output)
  --overwrite              Allow overwriting existing outputs
```

Chain selections accept comma-separated lists and ranges (e.g. `0,2,5-7`). Use
the interactive prompts when the `--receptor`/`--ligand` flags are omitted.
