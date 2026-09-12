# Using NebulaCell

This guide walks through everything the app can currently do, in the order you'll typically use it.

## 1. Starting the app

\`\`\`bash
cd Nebulacell/nebulacell
source .venv/bin/activate
python -m app.main
\`\`\`

Open `http://127.0.0.1:8050`. On first launch, a small demo dataset (Squidpy's Visium mouse brain, 2,800 cells) loads automatically — you can start exploring immediately without any data of your own.

## 2. Loading your own data

NebulaCell reads standard `.h5ad` (AnnData) files. **Requirement:** the file must have spatial coordinates stored in `adata.obsm['spatial']` — this is the standard convention for spatial transcriptomics data (Visium, MERFISH, Xenium, and most other spatial platforms all use it).

1. Type the path to your `.h5ad` file into the text box near the top (relative or absolute path — e.g. `data/anndata/my_sample.h5ad`).
2. Click **Load dataset**.
3. Wait for the label above to update with the cell count. Large datasets (see below) take longer.

If the file has no `spatial` key in `.obsm`, you'll get a clear error message rather than a crash.

### In-memory vs. "backed mode"

NebulaCell decides how to load a file based on its actual size (cells × genes), not just cell count. Small-to-medium datasets load fully into memory immediately. Very large datasets (currently: over ~50 million total values) are opened in a lighter "backed" mode first, then materialized into memory before any computation — so very large files take longer to load, and you'll see a `backed mode — large dataset` note in the label once loaded.

**Tested scale:** confirmed working from 2,800 cells up to 827,000+ cells (10x Xenium human lung).

## 3. Choosing a coloring mode

Three radio options control what the scatter plot's colors represent:

- **Leiden (expression clusters)** — uses expression-based clusters, if the loaded file already has a `leiden` column in `.obs`. If not present, this mode automatically falls back to spatial domain coloring instead of crashing.
- **HDBSCAN (spatial domains)** — colors cells by density-based spatial clustering, computed independently of gene expression, purely from physical coordinates. See section 5 for tuning this.
- **Gene expression** — colors every cell by the expression level of whichever gene is selected in the dropdown below.

## 4. Exploring gene expression

- The **Gene** dropdown lists genes available in the currently loaded dataset (up to 500, prioritizing highly-variable genes if that annotation exists in the file).
- Selecting a gene and switching to **Gene expression** coloring mode recolors every point by that gene's value.
- **Hovering any cell** always shows: the cell index, its top 3 most-expressed genes, and the currently selected gene's expression value — regardless of which coloring mode is active.

## 5. Tissue domain clustering (HDBSCAN)

Spatial domains are computed independently from expression clusters, based purely on physical density of points.

- Default `min_cluster_size` / `min_samples` values are **auto-scaled to your dataset's size** when it loads — larger datasets get larger defaults, since dense continuous scans (MERFISH, Xenium) behave very differently from Visium's sparse hex-spot grid.
- To tune manually: change the two number inputs and click **Recluster**. This re-runs HDBSCAN using only the already-loaded spatial coordinates — it's fast, and does not reload the dataset or touch the expression matrix.
- Feedback text below the button reports how many domains were found. Switch to **HDBSCAN (spatial domains)** coloring mode to see the result visually.
- **Tip:** if you see almost everything collapse into one giant domain plus noise, try lowering `min_cluster_size`. If you see too many tiny fragmented domains, raise it.

## 6. Histology image overlay

If the loaded file has a Visium-style histology image (`adata.uns['spatial']`), it renders as a background layer beneath the scatter points, correctly scaled so spots align with tissue structure.

- Toggle it on/off with the **Show histology image** checkbox — useful for comparing point patterns against raw tissue morphology, or for stripping the image away when it's visually busy.
- If the loaded dataset has no histology image (imaging-based assays like MERFISH/Xenium typically don't), the checkbox has no visible effect since there's nothing to show — this is expected behavior, not a bug.

## 7. Natural language query filtering

Type a query into the text box near the bottom and click **Run query** to highlight matching cells; non-matching cells dim to low opacity (rather than disappearing entirely, so you keep spatial context).

**Supported query patterns:**

| Example | Effect |
|---|---|
| `show cells with high Xkr4` | Highlights cells in the top 25% of expression for that gene |
| `show cells with low Xkr4` | Highlights cells in the bottom 25% |
| `highlight cells expressing high GeneName` | Same as above — `show`/`highlight` and `with`/`expressing` are interchangeable |
| `show cluster 3` | Highlights only cells in HDBSCAN spatial domain `3` |
| `show only domain 5` | Same pattern, works with `cluster` or `domain` |

Clicking **Reset** clears the filter and the query text box.

**Notes:**
- Gene names and cluster IDs are validated against the actual loaded dataset before filtering — if you typo a gene name or reference a cluster that doesn't exist, you'll get a specific error message (not a silent no-op or a crash).
- This is a deliberately narrow, rule-based (regex) parser — not a general LLM-backed natural language interface. It runs at zero cost and can't hallucinate an invalid filter, but it also won't understand queries outside its supported patterns. If a query doesn't match anything, you'll see a "couldn't understand that query" message with example phrasing.

## Troubleshooting

**"File not found"** — double check the path is relative to where you launched `python -m app.main` from (typically the `nebulacell/` inner directory), or use an absolute path.

**"...has no 'spatial' coordinates in .obsm"** — the file isn't spatial transcriptomics data, or its coordinates are stored under a non-standard key. NebulaCell currently only recognizes the standard `obsm['spatial']` convention.

**App feels slow after loading a huge dataset** — expected for very large files; clustering and per-cell gene-lookup precomputation run once at load time across the full dataset. See the "Known Limitations" section in the README regarding out-of-core computation.

**Gene dropdown shows unfamiliar identifiers (e.g. `ENSG00000001626`)** — some source files store Ensembl IDs instead of gene symbols. This reflects how the original file was prepared; NebulaCell doesn't currently remap IDs to symbols.
