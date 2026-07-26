"""
ATM Protein Conservation Explorer - Streamlit MVP

A Streamlit app for exploring evolutionary conservation of ATM protein
variants across vertebrate orthologs.
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


@st.cache_data
def load_csv(path: str) -> pd.DataFrame:
    """Load a CSV file into a DataFrame."""
    return pd.read_csv(path)


def check_columns(
    df: pd.DataFrame,
    required_columns: list,
    file_label: str,
) -> None:
    """Stop the app if required columns are missing."""
    missing = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing:
        st.error(
            f"'{file_label}' is missing required column(s): {missing}\n\n"
            f"Columns actually found: {list(df.columns)}"
        )
        st.stop()


def show_figure(column, path: str, caption: str) -> None:
    """Display a saved figure or show a warning if it is missing."""
    if os.path.exists(path):
        column.image(
            path,
            caption=caption,
            use_container_width=True,
        )
    else:
        column.warning(f"Missing figure: `{path}`")


# Page setup
st.set_page_config(
    page_title="ATM Protein Conservation Explorer",
    layout="wide",
)

st.title("🧬 ATM Protein Conservation Explorer")

st.write(
    "Explore how conserved each ATM missense variant position is across "
    "vertebrate orthologs, alongside ClinVar review status and gnomAD "
    "population-frequency evidence."
)


# Load data
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

check_columns(
    variant_df,
    REQUIRED_VARIANT_COLUMNS,
    VARIANT_CSV,
)

check_columns(
    position_df,
    REQUIRED_POSITION_COLUMNS,
    POSITION_CSV,
)


# Variant selector
st.sidebar.header("Select a variant")

variant_options = sorted(
    variant_df["protein_change"]
    .dropna()
    .unique()
)

selected_variant = st.sidebar.selectbox(
    "Protein change (e.g. S2R)",
    variant_options,
)

variant_row = (
    variant_df[
        variant_df["protein_change"] == selected_variant
    ]
    .iloc[0]
)

selected_position = int(variant_row["protein_position"])
selected_conservation = float(
    variant_row["conservation_percent"]
)
usable_species = int(
    variant_row["usable_species_count"]
)
alternate_species = int(
    variant_row["alternate_species_count"]
)


# Variant summary
st.subheader(f"Variant summary: {selected_variant}")

col1, col2, col3, col4 = st.columns(4)

col1.metric(
    "Conservation (%)",
    f"{selected_conservation:.1f}",
)

col2.metric(
    "Usable species",
    usable_species,
)

col3.metric(
    "Species with alternate residue",
    alternate_species,
)

col4.metric(
    "gnomAD status",
    str(variant_row["gnomad_match_status"]),
)


# ATM-wide conservation plot
st.subheader("ATM-wide conservation")

fig, ax = plt.subplots(figsize=(12, 3.5))

ax.plot(
    position_df["protein_position"],
    position_df["conservation_percent"],
    linewidth=0.7,
    color="#2b6cb0",
    label="Conservation by position",
)

position_match = position_df[
    position_df["protein_position"] == selected_position
]

if not position_match.empty:
    highlight_y = float(
        position_match["conservation_percent"].iloc[0]
    )
else:
    highlight_y = selected_conservation

ax.axvline(
    selected_position,
    color="#c53030",
    linestyle="--",
    linewidth=1,
    alpha=0.6,
)

ax.scatter(
    [selected_position],
    [highlight_y],
    color="#c53030",
    s=80,
    zorder=5,
    label=(
        f"Selected: {selected_variant} "
        f"(position {selected_position})"
    ),
)

ax.set_xlabel("Protein position")
ax.set_ylabel("Conservation (%)")
ax.set_ylim(0, 105)
ax.legend(loc="lower left", fontsize=8)

fig.tight_layout()
st.pyplot(fig)

plt.close(fig)


# Evidence interpretation
st.subheader("Evidence interpretation")

available_scores = pd.to_numeric(
    position_df["conservation_percent"],
    errors="coerce",
).dropna()

available_coverage = pd.to_numeric(
    variant_df["usable_species_count"],
    errors="coerce",
).dropna()

if not available_scores.empty:
    conservation_percentile = (
        available_scores.le(selected_conservation).mean()
        * 100
    )

    median_conservation = float(
        available_scores.median()
    )
else:
    conservation_percentile = float("nan")
    median_conservation = float("nan")

if not available_coverage.empty:
    coverage_percentile = (
        available_coverage.le(usable_species).mean()
        * 100
    )
else:
    coverage_percentile = float("nan")


if pd.isna(conservation_percentile):
    relative_interpretation = (
        "The position could not be ranked against the "
        "other ATM positions."
    )

elif conservation_percentile >= 90:
    relative_interpretation = (
        "This position is among the most highly conserved "
        "ATM positions included in the analysis."
    )

elif conservation_percentile >= 75:
    relative_interpretation = (
        "This position has relatively high conservation "
        "compared with other analyzed ATM positions."
    )

elif conservation_percentile >= 25:
    relative_interpretation = (
        "This position falls within the middle range of "
        "ATM conservation scores."
    )

else:
    relative_interpretation = (
        "This position is among the less conserved ATM "
        "positions included in the analysis."
    )


if (
    not pd.isna(conservation_percentile)
    and conservation_percentile >= 75
):
    evidence_takeaway = (
        "The result is consistent with evolutionary constraint "
        "at this position. A change at the residue may therefore "
        "warrant closer investigation alongside clinical, "
        "population, structural, and experimental evidence."
    )

elif (
    not pd.isna(conservation_percentile)
    and conservation_percentile < 25
):
    evidence_takeaway = (
        "The lower relative conservation provides weaker "
        "evolutionary evidence that this position is strongly "
        "constrained. This does not establish that the variant "
        "is benign."
    )

else:
    evidence_takeaway = (
        "The conservation result provides moderate or mixed "
        "evolutionary evidence. Additional evidence is needed "
        "before making any claim about variant significance."
    )


if pd.isna(conservation_percentile):
    percentile_text = "Unavailable"
else:
    percentile_text = f"{conservation_percentile:.0f}th percentile"

if pd.isna(coverage_percentile):
    coverage_text = (
        f"The score uses {usable_species:,} usable orthologs."
    )
else:
    coverage_text = (
        f"The score uses **{usable_species:,} usable orthologs**, "
        f"placing its alignment coverage around the "
        f"**{coverage_percentile:.0f}th percentile** among "
        f"the analyzed variants."
    )


st.info(
    f"""
### {selected_variant} at ATM residue {selected_position}

The human residue is conserved in
**{selected_conservation:.1f}%** of usable vertebrate orthologs.

Its conservation score falls around the
**{percentile_text}** among the ATM positions included in this dataset.

{relative_interpretation}

{coverage_text}

A total of **{alternate_species:,} orthologs** contain an alternate
residue at this aligned position.

**Evolutionary takeaway:** {evidence_takeaway}
"""
)

st.caption(
    "Conservation is one evidence type only. It cannot independently "
    "classify a variant as pathogenic, benign, or clinically actionable."
)


# Selected variant details
st.subheader("Selected variant details")

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

detail_columns = [
    column
    for column in preferred_columns
    if column in variant_df.columns
]

detail_table = variant_row[
    detail_columns
].to_frame(name="Value")

st.dataframe(
    detail_table,
    use_container_width=True,
)


# Reference figures
with st.expander("ATM-wide dataset overview"):
    st.write(
        "These figures summarize the complete conservation dataset "
        "rather than only the currently selected variant."
    )

    fig_col1, fig_col2, fig_col3 = st.columns(3)

    show_figure(
        fig_col1,
        FIGURE_RANKED,
        "Most and least conserved positions",
    )

    show_figure(
        fig_col2,
        FIGURE_COVERAGE,
        "Species coverage per position",
    )

    show_figure(
        fig_col3,
        FIGURE_ALT_RESIDUE,
        "Alternate residue distribution",
    )


# Methods, data sources, and limitations
with st.expander("Methods, data sources, and limitations"):
    st.markdown(
        """
### Conservation score formula

For each protein position:

```text
conservation_score   = matching_species_count / usable_species_count
conservation_percent = conservation_score × 100
````

`matching_species_count` is the number of aligned ortholog species whose
residue matches the human residue at that position.

`usable_species_count` excludes alignment gaps (`-`) and ambiguous
residues (`X`).

### Data sources

* Variant and clinical review data: ClinVar
* Population frequency data: gnomAD
* Ortholog protein sequences: NCBI Datasets
* Multiple-sequence alignment: MAFFT

### Interpretation

A highly conserved position may indicate that the residue has been
maintained throughout evolution and could be under stronger functional
constraint.

However, conservation alone does not demonstrate that a particular
amino-acid substitution damages protein function.

### Limitations

* Species sampling is uneven across positions, so the number of usable
  species and confidence in each score vary.
* All species are weighted equally, despite differences in evolutionary
  distance from humans.
* Alignment gaps and ambiguous residues are excluded rather than treated
  as evidence for or against conservation.
* The score measures exact residue identity rather than physicochemical
  similarity between amino acids.
* The tool reflects sequence conservation only.
* It does not directly include protein structure, functional domains,
  experimental assays, segregation evidence, or validated clinical
  outcomes.

### Disclaimer

This tool is for exploratory, research, and educational purposes only.

**It is not intended for clinical diagnosis, clinical variant
classification, or medical decision-making.**
"""
)


