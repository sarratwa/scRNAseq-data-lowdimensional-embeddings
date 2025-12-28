import cellxgene_census
import anndata as ad
import tiledbsoma as soma
from pathlib import Path
import time
import psutil
import os
from datetime import datetime

#-------------------------------------------------------------------------------
# Sample size = 3_000 
# Sampling summary:
# Cells × genes Final anndata shape: 3000 cells × 61497 genes
# Runtime: 14.65 minutes
# CPU RAM increase: 2551.8 MB
#-------------------------------------------------------------------------------
# Sample size = 10_000 
# Sampling summary
#  Cells × genes: Final anndata shape: 10000 cells × 61497 genes
# Tissue breakdown:
# tissue
# dorsolateral prefrontal cortex    2768
# telencephalon                      893
# cerebral cortex                    569
# cerebellum                         565
# middle temporal gyrus              455
# Name: count, dtype: int64
# Runtime: 40.75 minutes
# CPU RAM increase: 2836.3 MB
#-------------------------------------------------------------------------------

# memory usage
process = psutil.Process(os.getpid())

def cpu_mem_mb():
    return process.memory_info().rss / 1024**2

def log(msg):
    now = datetime.now().strftime("%H:%M:%S")
    print(f"[{now}] {msg}", flush=True)

OUTPUT_DIR = Path("data/raw")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)      

SAMPLE_SIZE = 10_000   # change to 3000 or 10_000
OUTPUT_FILE = OUTPUT_DIR / f"brain_{SAMPLE_SIZE}_sample.h5ad"

mem_start = cpu_mem_mb()
time_start = time.perf_counter()

print("Opening census...")
with cellxgene_census.open_soma(census_version="latest") as census:
    print("Connected")

    human = census["census_data"]["homo_sapiens"]

    print("Reading brain obs metadata...")
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

time_end = time.perf_counter()
mem_end = cpu_mem_mb()

log("Download complete")

print("\n--- Sampling summary ---")
print(f"\n Cells × genes: Final anndata shape: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Tissue breakdown:\n{adata.obs['tissue'].value_counts().head()}")
print(f"Runtime: {(time_end - time_start) / 60:.2f} minutes")
print(f"CPU RAM increase: {mem_end - mem_start:.1f} MB")

print("\nSaving...")
adata.write_h5ad(OUTPUT_FILE)
print(f"Saved to: {OUTPUT_FILE}")

'''
Chunk downloading for bigger samples
# config
SAMPLE_SIZE = 10_000          # change to 3000 if needed
CHUNK_SIZE = 2000             # 1000–2000 depending on performance
RANDOM_STATE = 42

OUTPUT_DIR = Path("data/raw")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / f"brain_{SAMPLE_SIZE}_sample.h5ad"

# memory tracking
process = psutil.Process(os.getpid())

def cpu_mem_mb():
    return process.memory_info().rss / 1024**2

def log(msg):
    now = datetime.now().strftime("%H:%M:%S")
    print(f"[{now}] {msg}", flush=True)

log(f"Sampling {SAMPLE_SIZE:,} brain cells (chunk size = {CHUNK_SIZE})")

mem_start = cpu_mem_mb()
time_start = time.perf_counter()


# open census
log("Opening CellxGene Census")
with cellxgene_census.open_soma(census_version="latest") as census:
    log("Connected")

    human = census["census_data"]["homo_sapiens"]

    # reading metadata
    log("Reading brain obs metadata")
    obs_reader = human.obs.read(
        value_filter="tissue_general == 'brain' and is_primary_data == True",
        column_names=[
            "soma_joinid", "cell_type", "tissue_general", "tissue",
            "donor_id", "assay", "disease", "sex"
        ]
    )

    obs_df = obs_reader.concat().to_pandas()
    log(f"Found {len(obs_df):,} candidate brain cells")

    # sampling cells
    n = min(SAMPLE_SIZE, len(obs_df))
    sampled = obs_df.sample(n=n, random_state=RANDOM_STATE)
    soma_ids = sampled["soma_joinid"].tolist()
    log(f"Selected {len(soma_ids):,} cells")

    # chunk download
    adatas = []
    n_chunks = (len(soma_ids) - 1) // CHUNK_SIZE + 1
    log(f"Downloading expression matrix in {n_chunks} chunks")

    for i in range(0, len(soma_ids), CHUNK_SIZE):
        chunk_ids = soma_ids[i:i + CHUNK_SIZE]
        chunk_idx = i // CHUNK_SIZE + 1

        log(f"→ Chunk {chunk_idx}/{n_chunks} ({len(chunk_ids)} cells)")

        with human.axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(coords=(chunk_ids,))
        ) as query:
            adata_chunk = query.to_anndata(
                X_name="raw",
                column_names={
                    "obs": [
                        "cell_type", "tissue_general", "tissue",
                        "donor_id", "assay", "disease", "sex"
                    ]
                }
            )

        adatas.append(adata_chunk)

# concat and save
log("Concatenating chunks")
adata = ad.concat(adatas, join="outer")

time_end = time.perf_counter()
mem_end = cpu_mem_mb()

log("Download complete")

print("\n--- Sampling summary ---")
print(f"Cells × genes: {adata.n_obs:,} × {adata.n_vars:,}")
print(f"Runtime: {(time_end - time_start) / 60:.2f} minutes")
print(f"CPU RAM increase: {mem_end - mem_start:.1f} MB")

log("Saving AnnData")
adata.write_h5ad(OUTPUT_FILE)
log(f"Saved to: {OUTPUT_FILE}")


'''
