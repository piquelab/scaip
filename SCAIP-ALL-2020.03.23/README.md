After align reads and deconvolute individuals, we would do the following analysis, including transforming data into seurat object, Batch clustering 
1. Read demux output, 1_demux2.R
2. Merge align reads from `/nfs/rprdata/julong/SCAIP/kallisto2/bus/, 2_merge_kb2.R
3. compare data of SCAIP1-5 and SCAIP1-6, 3_compareNewAndOld.R
4. Normalization data, correct batch effects and cluster analysis , 4_Harmony.R
5. Cell type annotation, 5_IdenCelltype_seurat.R
6. Differential expression analysis, 6_DEG.CelltypeNew.R
7. Gene ontology (GO) enrichment analysis for DEG, 7_GSE.ClusterProfiler.R
9. Calculating linear discriminant analysis (LDA), 9_RNA.dynamic2.R
10. Gene variability analysis,
1. 10_RNA.variance.R for calculating gene variability parameters and differential analysis;
2. 10_GSE.ClusterProfiler.R for enrichment analysis for DVG 
11. Example genes, 11_GENE.example.R
## 