import cellxgene_census
import anndata as ad
from pathlib import Path
import tiledbsoma as soma

# ---------- This file will have multiple commented sampling codes for 3000 cells and ----------
# ----------- Please chose the number of cells you will be working with ----------------

# -----------Data sampling of 10k cells -----------------
# --------- Runtime : ... minutes --------------------

# --------- Data sampling of 3000 cells -----------------
# --------- Runtime : 25 minutes --------------------
'''
OUTPUT_DIR = Path("data/processed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "brain_3000_sample.h5ad"
SAMPLE_SIZE = 3000  

print("Opening census...")
with cellxgene_census.open_soma(census_version="latest") as census:
    print("Connected")

    human = census["census_data"]["homo_sapiens"]

    print("Reading brain obs metadata (fast)...")
    obs_reader = human.obs.read(
        value_filter="tissue_general == 'brain' and is_primary_data == True",
        column_names=["soma_joinid", "cell_type", "tissue_general", "tissue", 
                      "donor_id", "assay", "disease", "sex"]
    )
    obs_df = obs_reader.concat().to_pandas()
    print(f"   → Found {len(obs_df):,} matching cells")

    # Safety limit
    n = min(SAMPLE_SIZE, len(obs_df))
    sampled = obs_df.sample(n=n, random_state=42)
    soma_ids = sampled["soma_joinid"].tolist()
    print(f"Selected {len(soma_ids):,} cells")

    print("⬇Downloading expression + obs using axis_query...")
    with human.axis_query(
        measurement_name="RNA",
        obs_query=soma.AxisQuery(coords=(soma_ids,))
    ) as query:
        adata = query.to_anndata(
            X_name="raw",
            column_names={
                "obs": ["cell_type", "tissue_general", "tissue",
                        "donor_id", "assay", "disease", "sex"]
            }
        )

print(f"\nFinal anndata shape: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Tissue breakdown:\n{adata.obs['tissue'].value_counts().head()}")
print("\nSaving...")
adata.write_h5ad(OUTPUT_FILE)
print(f"Saved to: {OUTPUT_FILE}")
'''