"""
ATM Protein Conservation Explorer - Streamlit MVP:

An mvp Streamlit app for exploring evolutionary
conservation of ATM protein variants across vertebrate orthologs.
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

# File paths
VARIANT_CSV = "data/processed/variant_scores.csv"
POSITION_CSV = "data/processed/position_scores.csv"

FIGURE_RANKED = "figures/ranked_conservation.png"
FIGURE_COVERAGE = "figures/species_coverage.png"
FIGURE_ALT_RESIDUE = "figures/alternate_residue_distribution.png"

REQUIRED_VARIANT_COLUMNS = [
    "protein_change",
    "protein_position",
    "conservation_percent",
    "usable_species_count",
    "alternate_species_count",
    "gnomad_match_status",
]

REQUIRED_POSITION_COLUMNS = [
    "protein_position",
    "conservation_percent",
]

# Data loading (cached so each CSV is only read once per session)
@st.cache_data
def load_csv(path: str) -> pd.DataFrame:
    """Load a CSV file into a DataFrame. Raises FileNotFoundError if missing."""
    return pd.read_csv(path)


def check_columns(df: pd.DataFrame, required_columns: list, file_label: str) -> None:
    """Stop the app with a clear message if required columns are missing."""
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        st.error(
            f"'{file_label}' is missing required column(s): {missing}\n\n"
            f"Columns actually found in the file: {list(df.columns)}"
        )
        st.stop()

# Page setup
st.set_page_config(page_title="ATM Protein Conservation Explorer", layout="wide")

st.title("🧬 ATM Protein Conservation Explorer")
st.write(
    "Explore how conserved each ATM missense variant position is across "
    "vertebrate orthologs, alongside ClinVar review status and gnomAD "
    "population-frequency evidence."
)

# Load data with error
try:
    variant_df = load_csv(VARIANT_CSV)
except FileNotFoundError:
    st.error(f"Could not find the variant data file: `{VARIANT_CSV}`")
    st.stop()

try:
    position_df = load_csv(POSITION_CSV)
except FileNotFoundError:
    st.error(f"Could not find the position data file: `{POSITION_CSV}`")
    st.stop()

check_columns(variant_df, REQUIRED_VARIANT_COLUMNS, VARIANT_CSV)
check_columns(position_df, REQUIRED_POSITION_COLUMNS, POSITION_CSV)

# Sidebar: variant selector
st.sidebar.header("Select a variant")

variant_options = sorted(variant_df["protein_change"].dropna().unique())

selected_variant = st.sidebar.selectbox("Protein change (e.g. S2R)", variant_options)

variant_row = variant_df[variant_df["protein_change"] == selected_variant].iloc[0]

# Metric cards for the selected variant
st.subheader(f"Variant summary: {selected_variant}")

col1, col2, col3, col4 = st.columns(4)
col1.metric("Conservation (%)", f"{variant_row['conservation_percent']:.1f}")
col2.metric("Usable species", int(variant_row["usable_species_count"]))
col3.metric("Species with alt. residue", int(variant_row["alternate_species_count"]))
col4.metric("gnomAD status", str(variant_row["gnomad_match_status"]))


# Dynamic ATM-wide conservation plot, highlighting the selected variant
st.subheader("ATM-wide conservation")

fig, ax = plt.subplots(figsize=(12, 3.5))

# The full conservation curve across every position that has a variant.
ax.plot(
    position_df["protein_position"],
    position_df["conservation_percent"],
    linewidth=0.7,
    color="#2b6cb0",
    label="Conservation by position",
)

selected_position = variant_row["protein_position"]

# Look up the position's own conservation value so the highlighted dot
# sits exactly on the curve (falls back to the variant's own value if the
# position isn't found in position_df for some reason).
position_match = position_df[position_df["protein_position"] == selected_position]
if not position_match.empty:
    highlight_y = position_match["conservation_percent"].iloc[0]
else:
    highlight_y = variant_row["conservation_percent"]

ax.axvline(selected_position, color="#c53030", linestyle="--", linewidth=1, alpha=0.6)
ax.scatter(
    [selected_position],
    [highlight_y],
    color="#c53030",
    s=80,
    zorder=5,
    label=f"Selected: {selected_variant} (pos {int(selected_position)})",
)

ax.set_xlabel("Protein position")
ax.set_ylabel("Conservation (%)")
ax.set_ylim(0, 105)
ax.legend(loc="lower left", fontsize=8)
fig.tight_layout()

st.pyplot(fig)

# Static figures generated during data processing
st.subheader("Reference figures")

fig_col1, fig_col2, fig_col3 = st.columns(3)


def show_figure(column, path, caption):
    """Display a saved PNG, or a clear warning if the file is missing."""
    if os.path.exists(path):
        column.image(path, caption=caption, use_container_width=True)
    else:
        column.warning(f"Missing figure: `{path}`")


show_figure(fig_col1, FIGURE_RANKED, "Most / least conserved positions")
show_figure(fig_col2, FIGURE_COVERAGE, "Species coverage per position")
show_figure(fig_col3, FIGURE_ALT_RESIDUE, "Alternate residue distribution")

# Selected variant details table
st.subheader("Selected variant details")

# Only include columns that actually exist in the loaded file
preferred_columns = [
    "gene",
    "protein_change",
    "hgvs_c",
    "hgvs_p",
    "transcript",
    "rs_id",
    "clinvar_variation",
    "review_status",
    "protein_position",
    "conservation_percent",
    "usable_species_count",
    "alternate_species_count",
    "gnomad_match_status",
    "gnomad_af",
]
detail_columns = [col for col in preferred_columns if col in variant_df.columns]

detail_table = variant_row[detail_columns].to_frame(name="Value")
st.dataframe(detail_table, use_container_width=True)

# Methods, dsta sources, and limitations
with st.expander("Methods, data sources, and limitations"):
    st.markdown(
        """
**Conservation score formula**

For each protein position:

```
conservation_score   = matching_species_count / usable_species_count
conservation_percent = conservation_score * 100
```

`matching_species_count` is the number of aligned ortholog species whose
residue matches the human residue at that position. `usable_species_count`
excludes alignment gaps ("-") and ambiguous residues ("X").

**Data sources**

- Variant and clinical review data: ClinVar
- Population frequency data: gnomAD
- Ortholog protein sequences: NCBI Datasets (vertebrate ATM orthologs)
- Multiple sequence alignment: MAFFT

**Limitations**

- Species sampling is uneven across positions, so the number of usable
  species (and therefore confidence in the score) varies.
- Alignment gaps and ambiguous residues are excluded rather than counted
  as evidence for or against conservation at that position.
- This tool reflects *sequence* conservation only — it does not account
  for protein structure, domain function, or experimental evidence.

**Disclaimer**

This tool is for exploratory and educational purposes only.
**It is not intended for clinical diagnosis or clinical decision-making.**
        """
    )
