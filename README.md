# ATM Protein Conservation Explorer

A small Streamlit app for exploring how evolutionarily conserved each ATM
missense variant position is across vertebrate orthologs, alongside ClinVar
review status and gnomAD population-frequency evidence.

## What it does

For a variant you pick from the sidebar, the app shows:

- Conservation percent, usable species count, alternate-residue species
  count, and gnomAD match status as metric cards
- A plot of conservation across the whole ATM protein, with the selected
  variant's position highlighted
- Three reference figures generated during data processing
- A details table for the selected variant
- An expander explaining the conservation score formula, data sources, and
  limitations

**This tool is for exploratory/educational use only — not for clinical
diagnosis or clinical decision-making.**

## How the data gets here

The app itself does no data processing. It just reads two CSVs that are
produced by an earlier notebook pipeline (not included in this folder):

1. Missense variants are pulled from ClinVar and cross-checked against
   gnomAD → `prepared_atm_missense_vus.csv`
2. Vertebrate ATM orthologs are downloaded and filtered →
   `atm_representative_proteins.fasta`
3. Orthologs are aligned with MAFFT →
   `atm_representative_proteins_aligned.fasta`
4. Each protein position's conservation score is computed from the
   alignment, then merged with the variant table → **`position_scores.csv`**
   and **`variant_scores.csv`** (these two are what the app reads)
5. A further-enriched version with extra tags (domain, conservation tier,
   evidence flag) is saved separately as `atm_variant_explorer_final.csv`
   — the app doesn't use this one, but it's there if you want to extend
   the app later.

Conservation score, per position:

```
conservation_score   = matching_species_count / usable_species_count
conservation_percent = conservation_score * 100
```

## Project layout the app expects

```
your-project-root/
├── app.py
├── requirements.txt
├── data/
│   └── processed/
│       ├── variant_scores.csv
│       └── position_scores.csv
└── figures/
    ├── ranked_conservation.png
    ├── species_coverage.png
    └── alternate_residue_distribution.png
```

## How to run it

```bash
pip install -r requirements.txt
streamlit run app.py
```

Run this from the project root (the folder containing `data/` and
`figures/`) so the relative paths in `app.py` resolve correctly.

## Notes for future me

- If a required CSV or column goes missing, the app fails loudly with a
  message listing what it expected vs. what it found — check
  `check_columns()` in `app.py` if that happens.
- The app is intentionally one file. If it grows (more pages, more data
  sources), that's the natural point to split it into modules — not
  before.