# Treehouse_OutlierRNASeq
## https://docs.google.com/presentation/d/1CQoIRCnTSLGydnxKiOXGFSrRBvqd_PYEHHUd6kU-9Jg/edit?usp=sharing
Undergraduate research project for the Treehouse Initiative by Liam McKay (undergrad) and Holly Beale (mentor)

Used cohort CKCC 146 Samples
All TH01 samples are Ribosomal RNA Depleted (RiboD)
TH02,3,4,6 are Polyadenylated RNA Selection (PolyA-S)

What we found:
There are high and low 95th percentiles -- Why?
Current findings:
- Grouping by sample center created more correlation between ckcc expression values and 
	- 95th percentile than pan-center
- RiboD samples have a lower 95th percentile than PolyA-S
- RiboD samples have a lower Number of Expressed Genes than PolyA-S



How to use this repository to generate my graphs:
- put ckcc_rsem_genes_results and comp4.3_tert8.ckcc.outlier_results folders in the same level as this README
- for all scatterplots (including grouping by TH0#) go to scatterplotsUMEND-NumOfExprGene-p95.R

- for all histograms use 
- single plot histogram: singleSamplePlot.R
- high and low 95th percentiles: GlobalPercentilePlots.R
- high and low 95th/75th percentiles: GlobalPercentilePlots.R

- all histograms of every sample: worstBestSamplesPlots.R
	- This also computes the "bump" value which is the highest point after the dip at 1.8 or 1.9 log2(TPM+1)
	- create a folder called BatchPlots

- high variance histograms of every sample: mostVariableGenes.R
	- create a folder called Batch-MostVariantGenesSorted-by-p95

- only top5 log2(TPM+1) genes histogram of every sample: top5thGenesPerSamples.R
	- create a folder called Batch-top5-sample-histograms

- expected count histogram batch: expectedCountScatter.R
	- create a folder called Batch-expectedCount

::==> create all of these folders in the same level as this file

- Boxplots for RiboD/PolyA-S 95th pctl and # of expressed gene comparison: boxplot.R
