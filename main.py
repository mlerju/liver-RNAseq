import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score

class RNASeqAnalyzer:
    def __init__(self, df, log_transform=True):
        """
        Parameters
        ----------
        df : pd.DataFrame
            Genes x Samples raw count matrix.
        log_transform : bool
            Apply log2(x+1) transformation.
        """
        self.df = df.copy()
        self.df_log = np.log2(df + 1) if log_transform else df.copy()
        self.labels = None
        self.marker_genes = {}
        self.enrichment_results = {}
        self.tsne_result = None

    # ----------------- Dimensionality Reduction -----------------
    def run_pca(self, n_components=50):
        self.df_t = self.df_log.T
        self.pca = PCA(n_components=n_components)
        self.pca_result = self.pca.fit_transform(self.df_t)
        return self.pca_result

    def run_tsne(self, perplexity=30, max_iter=1000):
        if not hasattr(self, 'pca_result') or self.pca_result.shape[1] < 50:
            self.run_pca(n_components=50)
        tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity, max_iter=max_iter)
        self.tsne_result = tsne.fit_transform(self.pca_result)
        return self.tsne_result

    def plot_tsne(self, labels=None, annotated_labels=None, figsize=(10,8), palette="tab10"):
        df_plot = pd.DataFrame({
            "tSNE1": self.tsne_result[:,0],
            "tSNE2": self.tsne_result[:,1],
            "Cluster": labels if labels is not None else self.labels
        })
        if annotated_labels is not None:
            df_plot["CellType"] = annotated_labels
            hue_col = "CellType"
        else:
            hue_col = "Cluster"
        plt.figure(figsize=figsize)
        sns.scatterplot(
            data=df_plot,
            x="tSNE1", y="tSNE2",
            hue=hue_col,
            palette=palette,
            s=50, alpha=0.8
        )
        plt.title("t-SNE Plot")
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show()

    # ----------------- Clustering -----------------
    def run_kmeans(self, n_clusters=3):
        self.kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        self.labels = self.kmeans.fit_predict(self.tsne_result)
        score = silhouette_score(self.tsne_result, self.labels)
        print(f"K-means silhouette score: {score:.3f}")
        return self.labels

    def run_dbscan(self, eps=5, min_samples=10):
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        labels = dbscan.fit_predict(self.tsne_result)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        print(f"DBSCAN found {n_clusters} clusters")
        if n_clusters > 1:
            score = silhouette_score(self.tsne_result, labels)
            print(f"DBSCAN silhouette score: {score:.3f}")
        return labels

    # ----------------- Marker Gene Detection -----------------
    def detect_marker_genes(self, top_n=50):
        if self.labels is None:
            raise ValueError("Clustering labels not found. Run clustering first.")
        self.marker_genes = {}
        for cluster in np.unique(self.labels):
            cluster_idx = np.where(self.labels == cluster)[0]
            other_idx = np.where(self.labels != cluster)[0]

            cluster_data = self.df_log.iloc[:, cluster_idx].values  # genes x cluster samples
            other_data = self.df_log.iloc[:, other_idx].values  # genes x other samples

            # Compute t-test for each gene individually
            t_stats = []
            p_vals = []
            for i in range(cluster_data.shape[0]):
                t, p = ttest_ind(cluster_data[i, :], other_data[i, :], equal_var=False)
                t_stats.append(t)
                p_vals.append(p)

            results = pd.DataFrame({
                "gene": self.df_log.index,
                "t_stat": t_stats,
                "p_val": p_vals,
                "mean_in_cluster": cluster_data.mean(axis=1),
                "mean_out_cluster": other_data.mean(axis=1)
            })
            results["log2FC"] = results["mean_in_cluster"] - results["mean_out_cluster"]
            results["adj_pval"] = multipletests(results["p_val"], method="fdr_bh")[1]  # FDR
            markers = results.sort_values(["adj_pval", "log2FC"], ascending=[True, False]).head(top_n)
            self.marker_genes[cluster] = markers
        return self.marker_genes

    # ----------------- Enrichment -----------------
    def run_enrichment(self, gene_sets=None, cutoff=0.05, top_n=200):
        gene_sets = gene_sets or ['KEGG_2021_Human', 'GO_Biological_Process_2021']
        self.enrichment_results = {}
        for cluster, markers in self.marker_genes.items():
            genes = markers['gene'].tolist()
            recognized_genes = [g for g in genes if not g.startswith('MT-') and len(g) <= 15]
            if len(recognized_genes) < 10:
                print(f"Cluster {cluster}: too few recognized genes for enrichment")
                continue
            enr = gp.enrichr(
                gene_list=recognized_genes,
                gene_sets=gene_sets,
                organism='Human',
                outdir=None,
                cutoff=cutoff
            )
            if not enr.results.empty:
                enr.results.rename(columns=str.title, inplace=True)
                self.enrichment_results[cluster] = enr.results
                print(f"Cluster {cluster} enrichment done, top terms:")
                print(enr.results[['Term','Adjusted P-Value']].head(5))
        return self.enrichment_results

    # ----------------- Enrichment Plotting -----------------
    def plot_enrichment(self, top_n=10, kind="bar", multi_cluster=False):
        """
        Plot enrichment results. Handles single or multi-cluster plotting.
        """
        def _standardize_columns(df):
            df = df.copy()
            df.columns = [c.strip().title().replace(" ", "_").replace("-", "_") for c in df.columns]
            adj_p_col = next((c for c in df.columns if "Adjusted" in c and "P" in c), None)
            score_col = next((c for c in df.columns if "Combined" in c and "Score" in c), adj_p_col)
            return df, adj_p_col, score_col

        if multi_cluster:
            all_data = []
            for cluster_id, df in self.enrichment_results.items():
                if df is None or df.empty:
                    continue
                df, adj_p_col, score_col = _standardize_columns(df)
                if not adj_p_col:
                    continue
                top_df = df.nsmallest(top_n, adj_p_col).copy()
                top_df["Cluster"] = cluster_id
                top_df["-log10P"] = -np.log10(top_df[adj_p_col])
                top_df.rename(columns={score_col: "CombinedScore"}, inplace=True)
                all_data.append(top_df)
            if not all_data:
                print("No enrichment data to plot.")
                return
            plot_df = pd.concat(all_data)
            plt.figure(figsize=(10,6))
            sns.scatterplot(
                data=plot_df,
                x="Cluster",
                y="Term",
                size="CombinedScore",
                hue="-log10P",
                sizes=(50, 300),
                palette="viridis",
                alpha=0.8
            )
            plt.title(f"Top {top_n} Enriched Pathways Across Clusters")
            plt.xlabel("Cluster")
            plt.ylabel("Pathway")
            plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", title="Legend")
            plt.tight_layout()
            plt.show()
            return

        # Single-cluster plotting
        for cluster_id, df in self.enrichment_results.items():
            if df is None or df.empty:
                continue
            df, adj_p_col, score_col = _standardize_columns(df)
            if not adj_p_col:
                continue
            top_df = df.nsmallest(top_n, adj_p_col)

            if kind == "bar":
                plt.figure(figsize=(8,6))
                sns.barplot(
                    data=top_df,
                    y="Term",
                    x=score_col,
                    palette="viridis"
                )
                plt.title(f"Top {top_n} Enriched Pathways - Cluster {cluster_id}")
                plt.xlabel("Combined Score")
                plt.ylabel("")
                plt.tight_layout()
                plt.show()
            elif kind == "dot":
                plt.figure(figsize=(8,6))
                sns.scatterplot(
                    data=top_df,
                    x=-np.log10(top_df[adj_p_col]),
                    y="Term",
                    size="Overlap" if "Overlap" in top_df.columns else score_col,
                    hue=score_col,
                    palette="viridis",
                    sizes=(50,300),
                    alpha=0.7
                )
                plt.title(f"Top {top_n} Enriched Pathways - Cluster {cluster_id}")
                plt.xlabel("-log10(Adjusted P-value)")
                plt.ylabel("")
                plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
                plt.tight_layout()
                plt.show()


# ----------------- Load data -----------------
df = pd.read_csv("human_liver.tsv", sep="\t", index_col=0)
print(df.shape)
print(df.head())

# ----------------- Initialize analyzer -----------------
analyzer = RNASeqAnalyzer(df)

# ----------------- Dimensionality Reduction -----------------
analyzer.run_pca(n_components=50)
analyzer.run_tsne(perplexity=30, max_iter=1000)

# ----------------- Clustering -----------------
analyzer.run_kmeans(n_clusters=3)

# ----------------- Marker Gene Detection -----------------
analyzer.detect_marker_genes(top_n=50)

# ----------------- Gene Set Enrichment -----------------
gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2021']
analyzer.run_enrichment(gene_sets=gene_sets, cutoff=0.05)

# ----------------- t-SNE Plot -----------------
# Optional: provide cluster annotation
cluster_annotations = {
    0: "Hepatocytes",
    1: "Kupffer cells",
    2: "Endothelial cells",
    3: "Stellate cells",
    4: "Cholangiocytes",
    5: "Immune cells"
}
labels_annotated = [cluster_annotations.get(x, "Unknown") for x in analyzer.labels]
analyzer.plot_tsne(annotated_labels=labels_annotated)

# ----------------- Enrichment Plots -----------------
analyzer.plot_enrichment(top_n=10, kind="bar")             # per-cluster barplots
analyzer.plot_enrichment(top_n=5, multi_cluster=True)     # multi-cluster comparison dotplot
