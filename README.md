# Liver RNA-Seq Clustering and Pathway Analysis

This project explores publicly available RNA-seq data of **human liver samples** (903 samples × 35,238 genes)  
to identify transcriptomic patterns and enriched biological pathways.

### 🔬 Methods
1. **Data Preprocessing**
   - Log2 transformation and scaling
   - PCA and t-SNE for dimensionality reduction

2. **Unsupervised Clustering**
   - K-means clustering (`n_clusters=3`, silhouette score = 0.57)
   - DBSCAN for density-based structure detection (silhouette score = 0.54)

3. **Cluster Characterization**
   - Top 200 marker genes extracted per cluster
   - Gene set enrichment using KEGG pathways (via `gseapy.enrichr`)

4. **Visualization**
   - PCA and t-SNE scatterplots
   - Cluster-specific enrichment barplots
   - Multi-cluster summary dotplot

### 📊 Key Results
- **Cluster 0**: Enriched for *complement and coagulation cascades*, *drug metabolism (CYP450)*.
- **Cluster 1**: Enriched for *ribosome*, *glycolysis/gluconeogenesis*, *HIF-1 signaling*.
- **Cluster 2**: Enriched for *ribosome*, *Salmonella/E. coli infection pathways*.

These findings highlight distinct transcriptional programs in human liver, capturing immune, metabolic, and stress-response functions.

### ⚙️ Technologies
- **Python**: pandas, scikit-learn, matplotlib, seaborn
- **Bioinformatics**: gseapy (Enrichr)
- **Clustering**: K-means, DBSCAN
- **Dimensionality Reduction**: PCA, t-SNE

### 🚀 Future Directions
- Add batch-effect correction and normalization
- Compare clustering on TPM vs raw counts
- Apply pathway network visualization (EnrichmentMap, Cytoscape)

---

## 💡 Takeaways
This project demonstrates:
- Handling **high-dimensional transcriptomics data**
- Applying **unsupervised learning** for biological discovery
- Performing **pathway enrichment** to interpret clusters
