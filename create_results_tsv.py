import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# Constants for ICTV taxonomy columns
TAXONOMY_COLUMNS = [
    "Realm",
    "Subrealm",
    "Kingdom",
    "Subkingdom",
    "Phylum",
    "Subphylum",
    "Class",
    "Subclass",
    "Order",
    "Suborder",
    "Family",
    "Subfamily",
    "Genus",
    "Subgenus",
    "Species",
]


# Apply the taxonomic classification based on the `gani` value
def assign_taxonomy(row):
    if row["gani"] >= 0.95:
        # Assign the gANI value to Species_score and keep Species
        row["Species_score"] = row["gani"]
    else:
        # If gani < 0.95, blank Species and find
        # the last non-empty taxonomic level
        row["Species"] = np.nan
        row['Subgenus'] = np.nan
        # Start from 'Genus' backwards
        taxonomic_columns = TAXONOMY_COLUMNS[-3::-1]
        score_columns = [f"{col}_score" for col in taxonomic_columns]
        for taxon, score_col in zip(taxonomic_columns, score_columns):
            if pd.notna(row[taxon]):
                row[score_col] = row["gani"]
                break  # Stop once the first non-empty taxon is found
    return row


# Helper function to display section headers
def print_section(step_num, title):
    print(f'\n[Step {step_num}] {title}\n{"-" * 50}')


# Create the argument parser
parser = argparse.ArgumentParser(
    description="Classify ICTV challenge contigs into taxonomy using "
    "global ANI (gANI) from Vclust output."
)
# Input ANI file
parser.add_argument(
    "--input-ani",
    type=Path,
    required=True,
    help="Path to the input ANI file",
)
# Input ICTV taxonomy file
parser.add_argument(
    "--input-taxonomy",
    type=Path,
    required=True,
    help="Path to the reference ICTV taxonomy file",
)
# Output file for results
parser.add_argument(
    "--output-tsv",
    type=str,
    default="vclust_results.tsv",
    help="Path to the output TSV file (default: %(default)s).",
)

# Parse arguments
args = parser.parse_args()


# 1. Load ICTV taxonomy
print_section(1, "Loading Reference ICTV Taxonomy")
print(f"- Input file: {args.input_taxonomy}")
taxonomy_data = pd.read_csv(
    args.input_taxonomy, sep="\t", usecols=["accession"] + TAXONOMY_COLUMNS
)
print(f"  * Genome accessions loaded: {len(taxonomy_data):,}")

# 2. Load Vclust ANI data
print_section(2, "Loading Vclust Output")
print(f"- Input file: {args.input_ani}")
ani_data = pd.read_csv(
    args.input_ani, sep="\t", usecols=["query", "reference", "gani"]
)
print(f"  * Total genome pairs loaded: {len(ani_data):,}")

# 3. Filter Vclust ANI data
print_section(3, "Filtering Genome Pairs")

# 3.1 Select only challenge-to-database hits
print("1. Challenge-to-database hits:")
ani_data = ani_data[
    ani_data["query"].str.startswith("ICTV")
    & ~ani_data["reference"].str.startswith("ICTV")
]
print(f"   * Genome pairs retained: {len(ani_data):,}\n")

# 3.2 Apply gANI filter (>= 70%)
print("2. Applying gANI filter (>= 70%):")
ani_data = ani_data[ani_data["gani"] >= 0.7]
print(f"   * Genome pairs retained: {len(ani_data):,}\n")

# 3.3 Best database hit per challenge contig
print("3. Best database hit per challenge contig:")
ani_data = ani_data.loc[ani_data.groupby("query")["gani"].idxmax()]
print(f"   * Genome pairs retained: {len(ani_data):,}")

# 4. Add taxonomy lineage to the data
print_section(4, "Assigning contigs to ICTV Taxonomy")
out_data = ani_data.merge(
    taxonomy_data, left_on="reference", right_on="accession", how="left"
)

# 4.1 Prepare output columns
output_columns = ["query", "reference", "gani"] + [
    item for col in TAXONOMY_COLUMNS for item in (col, f"{col}_score")
]
out_data = out_data.reindex(columns=output_columns)

# 4.2 Update scores for species and genus based on gANI values
out_data = out_data.apply(assign_taxonomy, axis=1)

# 4.3 Clean up dataframe (drop unnecessary columns)
out_data.drop(columns=["reference", "gani"], inplace=True)
out_data.rename(columns={"query": "SequenceID"}, inplace=True)

# 4.4 Summary of predictions
print(f"- Total Contigs: {len(out_data):,}")
for tax_rank in TAXONOMY_COLUMNS[::-1]:
    count = out_data[f"{tax_rank}_score"].count()
    if count:
        print(f"   * {tax_rank:<10}: {count:,}")

# 5. Save results to CSV
print_section(5, "Saving Results")
print(f"- Output file: {args.output_tsv}")
# Sort results by SequenceID
out_data['tmp'] = out_data['SequenceID'].str.extract(r'_(\d+)$').astype(int)
out_data = out_data.sort_values(by='tmp').drop(columns=['tmp'])
# Replace NaN with 'NA'
out_data = out_data.astype("object").fillna("NA")
out_data.to_csv(args.output_tsv, sep='\t', index=False)
print("- Status     : Done! Results saved.")
# print()