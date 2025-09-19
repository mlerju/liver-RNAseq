# Liver RNA-seq Clustering & Marker Gene Analysis

This repository contains a **Python pipeline** for clustering human liver RNA-seq samples, identifying cluster-specific marker genes, and performing pathway enrichment analysis. The project includes dimensionality reduction (PCA, t-SNE), clustering (K-means, DBSCAN), differential expression, and gene set enrichment with visualization.

---

## Features

- **Preprocessing & Normalization**
  - Log-transformation of raw counts
  - Filtering lowly expressed genes

- **Dimensionality Reduction**
  - PCA for variance exploration
  - t-SNE for visualization

- **Clustering**
  - K-means
  - DBSCAN
  - Silhouette score evaluation

- **Marker Gene Detection**
  - Gene-wise t-test (cluster vs rest)
  - Log2 fold-change and FDR-adjusted p-values
  - Top N marker genes per cluster

- **Gene Set Enrichment**
  - Uses [gseapy](https://github.com/zqfang/GSEApy) `enrichr` for KEGG and GO pathways
  - Filtering mitochondrial/unrecognized genes
  - Barplots and dotplots of enriched pathways

- **Visualization**
  - PCA and t-SNE scatterplots
  - Heatmaps of top marker genes
  - Enrichment barplots and multi-cluster dotplots

---

## Installation

1. Clone the repository:

git clone https://github.com/yourusername/liver-rnaseq-clustering.git
cd liver-rnaseq-clustering

    Create a Python virtual environment and activate it:

python -m venv .venv
source .venv/bin/activate    # Linux/macOS
.venv\Scripts\activate       # Windows

    Install dependencies:

pip install -r requirements.txt

Dependencies include:

    pandas

    numpy

    matplotlib

    seaborn

    scipy

    scikit-learn

    gseapy

Usage

    Prepare your RNA-seq dataset as a TSV file with genes as rows and samples as columns:

gene    sample1    sample2    sample3 ...
GENE1   100        50         200
GENE2   0          10         5
...

    Edit main.py to point to your dataset:

df = pd.read_csv("human_liver.tsv", sep="\t", index_col=0)

    Run the analysis:

python main.py

The pipeline will:

    Normalize the data

    Run PCA and t-SNE

    Perform K-means and DBSCAN clustering

    Detect cluster-specific marker genes

    Perform pathway enrichment

    Plot heatmaps, barplots, and multi-cluster enrichment visualizations

Project Structure

liver-rnaseq-clustering/
│
├── main.py                 # Main script for analysis
├── liver_rnaseq_analyzer.py  # Class-based modular pipeline
├── human_liver.tsv         # Example RNA-seq dataset
├── requirements.txt        # Python dependencies
└── README.md

Configuration

You can modify the following parameters in main.py:

    Clustering

kmeans_clusters = 3
dbscan_eps = 5
dbscan_min_samples = 10

Marker genes

top_n_markers = 50

Enrichment

gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2021']
enrichment_cutoff = 0.05

t-SNE

    perplexity = 30
    n_iter = 1000

Visualization

    PCA and t-SNE plots: clusters and sample distribution

    Heatmaps: top marker genes per cluster

    Barplots: top enriched pathways per cluster

    Dotplots: multi-cluster enrichment comparison

All plots are displayed interactively and can be saved by adding plt.savefig("filename.png") before plt.show().
Notes

    Mitochondrial (MT-) genes are filtered during enrichment to avoid bias.

    Adjust top_n_markers and cutoff to control marker gene selection and enrichment sensitivity.

    Ensure gseapy is installed and has internet access to query Enrichr.
