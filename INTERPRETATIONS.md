# Discussion of Data and Biological Interpretaions

> Agarwal, D.; Kuhns, R.; Dimitriou, C.N.; Barlow, E.; Wahlin, K.J.; Enke, R.A. Bulk RNA Sequencing Analysis of Developing Human Induced Pluripotent Cell-Derived Retinal Organoids. *Sci Data *2022*, 9*, 759, https://doi.org/10.1038/s41597-022-01853-x.

----

## Understanding the Experiment
Agarwal et. al. aimed to develop viable hiPSC derived retinal organoids with fluorescent reporters (SIX6 and CRX) to be later utilized to understand the development of the human eye.

They ran a time-series experiment with 7 time points, each point representing a different stage of retinal development. They had a total of 28 samples sequenced.

| Timepoint | Stage of Retinal Development | # of Replicates|
| --- | ---| ---|
| D0 | Control; hiPSCs| 3 |
| D10| Early Retinal Induction| 4 |
| D25| Late Retinal Induction| 4 |
| D65| PR Generation| 4 |
| D100| PR Generation with early RPE formation| 4 |
| D180| Mature Organoids with RPE presence| 6|
| D280| Late Stage  Organoids with mature RPE  | 3|

## Differential Expression Analysis
Ideally, pairwise contrasts would be drawn to indicate differentially expressed (DE) genes. While edgeR shines with small replicate numbers modelling batch effects, DESeq2 was chosen for analysis because of reproducibility and its fold change shrinkage.

The [PCA](https://github.com/jroverflow/repro-bulkrnaseq-agarwal2022/blob/main/fig3/panelA-PCAPlot.png) plot shows clear and separate clustering amongst timepoints D0-D100 samples. We can assume the PC1 variable is likely development separator. This is supported by how D180 and D280 samples are clustered around the same area. Further genomic differences between D180 and D280 organoids may likely be more thoroughly displayed in the PC3/4 dimensions.

This is further supported by the [sample-to-sample distance heatmap](https://github.com/jroverflow/repro-bulkrnaseq-agarwal2022/blob/main/fig3/panelB-stsHeatmap.png).There is clear clutering with the D0, D10, and D25 samples. The heatmap also clarifies the clutering of D180 and D280 samples in the PCA plot, as sample D280-A2 is shown to have the shortest distance to the D180-6 sample.

The contrast of D25 and D65 was chosen because it represents the commitment of the retinal progenitor into becoming a photoreceptor. Errors in early PR generation at D65 will likely not be fixed as the PR matures into an organoid (D180/D280). Its related [volcano plot](https://github.com/jroverflow/repro-bulkrnaseq-agarwal2022/blob/main/fig3/panelD-volcanoPlot.png) indicates upregulation of key retinal genetic markers at this critical point. This is also highlighted by the authors:

> This analysis demonstrates differential expression of early retinal progenitor transcripts including LIN28A, NES, MKI67, GMNN and GNL3 at the day 25 time point and maturing retina genes POU4F2, VSX2, CRX, ATOH7, ARR3 and LHX4 at day 65. (Agarwal et al. 2022).

It is important to note though that in the author's original volcano plot, it shows differential expression, in addition downregulation of maturing retinal markers at the D65, via the signage of log2foldchange. This is likely due to an incorrect order of contrast in R.

```R
# potentially this
results <- results(dds, contrast = c("Timepoint", "D25", "D65"))
# instead of this
results <- results(dds, contrast = c("Timepoint", "D65", "D25"))
```

## Pathway Enrichment Analysis
Overrepresentation analysis (ORA) of DE genes between the same contrast (D25 v D65) finds several promising signs towards viable retinal organoids. 

[BP of GO](https://github.com/jroverflow/repro-bulkrnaseq-agarwal2022/blob/main/additionalPlots/goEnrichmentPlot.png) enrichment analysis supports this claim. Explicitly, within the top 15 GO pathways, *sensory system development* (GO:0048880), *visual system development* (GO:0150063), and *eye development* (GO:0001654) indicate a significant (padj = 1e-09) portion of DE genes are expressed in those pathways. In addition, given that this contrast in particular is at the turning point of retinal progenitors into photoreceptors, the top two GO pathways: *regulation of trans-synaptic signaling* (GO:0099177) and *modulation of chemical synaptic transmission* (GO:0050804) are especially promising.

[KEGG enrichment](https://github.com/jroverflow/repro-bulkrnaseq-agarwal2022/blob/main/additionalPlots/keggEnrichmentPlot.png) analysis builds on the story from GO. Two particular pathways should be highlighted in relation to retinal development: axon guidance (GO:0007411) and neurotransmitter *glutamatergic synapse* (GO:0098978) are significantly overrepresented.

[GSEA](https://github.com/jroverflow/repro-bulkrnaseq-agarwal2022/blob/main/additionalPlots/gsea_top20pathways.png) plot rounds out story. Now considering all genes, a majority of the top 10 most upregulated GO pathways relate to retinal development and the visual system. All three levels of enrichment analysis indicate a positive look for the author's retinal organoids.
