# Methylation-Analysis
Reduced representation bisulfite sequencing (RRBS) was used to profile DNA methylation in human astrocytes in the context of ATP13A2 loss-of-function (LOF). Astrocytes were derived from isogenic human cell lines representing three experimental conditions: (i) wild-type (WT) isogenic controls, (ii) ATP13A2 c.1306 loss-of-function mutant cells (c.1306), and (iii) ATP13A2 c.1306 cells treated with CGP (c.1306+CGP). Each condition was represented by approximately three independent biological replicates. The primary goal of the computational analysis was to determine whether ATP13A2 LOF is associated with global or regional alterations in CpG methylation, in light of prior evidence that ATP13A2 dysfunction upregulates the enzyme AMD1, shifting S-adenosylmethionine (SAM) toward decarboxylated SAM (dcSAM) and potentially impacting methyl donor availability for DNA and histone methylation.




# RRBS methylation analysis (human astrocytes; ATP13A2 LOF ± CGP)

## Overview
This repository contains a generalized workflow for RRBS methylation analysis:
- Trimming + QC (CutAdapt, FastQC, MultiQC)
- Alignment + methylation calling (Bismark, Samtools, MethylDackel)
- Downstream methylation analysis in R (methylKit):
  - QC summaries and distributions
  - CpG-level and tiled (regional) methylation
  - Differential methylation with logistic regression + overdispersion correction
  - Export of DMR coordinates for UCSC Table Browser gene annotation
  - Optional integration with gene expression and enrichment (clusterProfiler)

## Experimental design
Conditions (approx. 3 biological replicates each):
- WT (isogenic control)
- ATP13A2 c.1306 LOF (c1306)
- c1306 + CGP (c1306_CGP)

## Inputs
### Raw sequencing data
- `*.fastq.gz` RRBS reads (paired or single-end depending on library)

### Methylation call files
Downstream R code expects CpG methylation calls in a methylKit-supported format, e.g.:
- Bismark coverage files (`*.bismark.cov` / `*.cov`, often gzipped)

Place these under `data/` and fill out:
- `metadata/samplesheet.csv`

## Preprocessing pipeline (command line outline)
> The R script starts AFTER methylation call files exist.

1) Trim adapters/primers/polyA with CutAdapt  
2) FastQC after trimming; summarize with MultiQC  
3) Align with Bismark to hg38; review Bismark reports  
4) Samtools alignment QC; confirm spike-in representation if present  
5) MethylDackel conversion efficiency (non-CpG + spike-ins) and methylation calls  
6) Collect methylation call files for each sample (coverage/bedGraph)

## Running the R analysis
From project root:

```bash
Rscript analysis_rrbs_methylkit.R



## Additional Information about the analysis

RRBS preprocessing, alignment, and raw data quality control
Raw RRBS reads in FASTQ format were processed using a standard bisulfite sequencing pipeline. Adapter sequences, primers, poly-A tails, and other unwanted bases at read ends were trimmed using CutAdapt. Trimming ensured removal of sequencing adapters and low-quality bases that could interfere with bisulfite alignment and methylation calling.
Post-trimming quality control was performed on each sample using FASTQC. FASTQC reports were inspected to evaluate: (i) per-base sequence quality, (ii) per-sequence quality score distributions, (iii) per-base GC content, (iv) sequence length distributions, and (v) the presence of overrepresented sequences or residual adapter contamination. Samples with clear evidence of poor base quality or persistent adapter contamination were flagged for closer inspection before downstream analysis.
Trimmed reads were aligned to the appropriate human reference genome using Bismark, a bisulfite-aware short-read aligner. Bismark performs directional alignment of bisulfite-converted reads and simultaneously records cytosine methylation status at CpG sites. For each sample, Bismark alignment reports were used to quantify (i) overall alignment efficiency, (ii) proportion of uniquely mapped reads, (iii) strand bias, and (iv) genome-wide coverage over CpG sites.
Alignment quality and coverage were further evaluated using Samtools. Samtools was used to inspect BAM files, calculate basic alignment statistics, and verify proper incorporation and coverage of spike-in control sequences. These steps ensured that coverage was sufficient and that no sample exhibited anomalous mapping patterns that would confound methylation estimates.
Bisulfite conversion efficiency was quantified using MethylDackel, which calculates methylation ratios at spike-in cytosines and non-CpG cytosines in the genome. Conversion rates were assessed at both non-CpG cytosines and spike-ins, and samples consistently exhibited >99% bisulfite conversion efficiency, supporting the reliability of observed CpG methylation levels.
Methylation calling, file generation, and import into R/methylKit
CpG-resolved methylation calls were generated during the Bismark/MethylDackel processing steps and summarized in tabular methylation call files. These files were made available via the RRBS MultiQC report and downloaded in bedGraph or BAM-derived formats. For very large files, downloads were carried out sequentially using standard command-line tools (e.g. curl) under Linux or Windows PowerShell to ensure robust transfer.
All downstream analyses were performed in R using the methylKit package. Methylation call files were imported with methylKit::methRead(), specifying sample identifiers, treatment groups, and genome assembly. Each imported sample was represented as a MethylRaw object containing CpG-level information, including chromosome, base position, strand, read coverage, and counts of methylated and unmethylated cytosines. Individual MethylRaw objects were then combined into a MethylRawList object to facilitate joint processing.
A typical CpG record contained the following fields: chromosome (chr), genomic coordinate (start/end), strand, total coverage (number of reads overlapping the cytosine), number of methylated cytosines, and number of unmethylated cytosines. This structure provided base-resolution methylation data suitable for both global and regional analyses.
Quality assessment and descriptive methylation statistics
A series of descriptive and QC analyses were performed to rigorously evaluate methylation data quality prior to inferential modeling. General sequencing and alignment statistics (e.g. total read counts, proportion of uniquely aligned reads, bisulfite conversion rates from spike-ins, total number of CpGs covered, and average CpG coverage) were compiled for each sample. Typical libraries yielded coverage of ~10–15 million unique CpG sites per sample, with mean CpG depth around 11–15×, indicating adequate depth for quantitative methylation analysis.
Within methylKit, the functions getMethylationStats() and getCoverageStats() were used to summarize and visualize CpG methylation and coverage distributions. getMethylationStats() generated histograms and summary statistics for the distribution of methylation percentages across CpG sites in each sample (e.g. minimum, quartiles, median, mean, maximum, and upper percentiles up to the 99th–100th percentile). getCoverageStats() provided analogous summaries for read coverage per CpG. Violin plots and histograms of these distributions were used to visually confirm that (i) methylation values spanned the expected 0–100% range without extreme artifacts, and (ii) coverage distributions were reasonably balanced across samples, with no evidence of systematic under- or over-coverage that might bias differential analyses. Together, these QC metrics supported inclusion of all samples in subsequent analyses.
Construction of a unified methylation object and global methylation quantification
To perform comparative analyses between conditions, per-sample MethylRaw objects were merged into a joint MethylBaseobject using the unite() function in methylKit. This procedure identified CpG loci that were commonly covered across samples (subject to user-defined minimum coverage and sample presence thresholds) and produced a single table in which each row corresponded to a CpG site and each column contained methylated and unmethylated read counts for each sample.
Global CpG methylation levels were estimated for each sample using the percMethylation() function, which calculates methylation percentages at each CpG site within the MethylBase object. Sample-wise global methylation values were then obtained by taking the mean methylation percentage across CpGs using colMeans(). These global averages were used to compare overall CpG methylation levels between ATP13A2 c.1306 astrocytes and WT controls, as well as to examine the effect of CGP treatment within the c.1306 background. The global comparisons were treated as descriptive summaries to contextualize more detailed regional analyses.
Regional methylation analysis using genomic tiling
To analyze spatial patterns of methylation beyond single CpGs, the genome was partitioned into fixed-width contiguous tiles using methylKit’s tileMethylCounts() function. This function aggregates methylation counts (methylated and unmethylated reads) for all CpGs falling within each genomic window, producing a new MethylRawList object in which each entry corresponds to a genomic tile rather than a single nucleotide. Tiling increases statistical power by integrating information across nearby CpGs and enables detection of differentially methylated regions (DMRs) that span multiple CpGs.
Tiles were generated across the genome for all samples, and the resulting tiled methylation data were converted into a MethylBase object for comparative analysis. For each tile (i.e. candidate DMR), the object contained methylated and unmethylated counts per sample, along with coverage and percent methylation values.
Statistical modeling of differential methylation and overdispersion correction
Differential methylation between experimental conditions was assessed using the logistic regression framework implemented in methylKit’s calculateDiffMeth() function. For each CpG site or genomic tile, let nini​ denote the coverage (number of reads overlapping that feature in sample ii) and πiπi​ the underlying methylation proportion. The number of methylated reads methimethi​ was modeled as a binomial random variable,
methi∼Binomial(ni,πi),methi​∼Binomial(ni​,πi​),
with expected mean μi=niπiμi​=ni​πi​ and variance Var(methi)=niπi(1−πi)Var(methi​)=ni​πi​(1−πi​) under the ideal binomial assumption. In practice, biological and technical variation often leads to overdispersion, whereby the observed variance of methylated counts exceeds niπi(1−πi)ni​πi​(1−πi​).
To account for this extra-binomial variability, methylKit estimates a scaling parameter ϕϕ that inflates the nominal variance. Overdispersion is quantified by comparing the residual deviance χ2χ2 of the fitted model to its degrees of freedom, and the scaling factor is computed as
ϕ=χ2/(N−P),ϕ=χ2/(N−P),
where NN is the number of observations and PP is the number of model parameters. This dispersion factor is then used to correct variance estimates, yielding more conservative and robust standard errors and p-values for the regression coefficients.
For inferential testing, methylation proportions were modeled on the logit scale. For a given site or region, the full logistic regression model included an indicator variable for treatment or genotype,
log⁡(πi1−πi)=β0+β1Treatmenti,log(1−πi​πi​​)=β0​+β1​Treatmenti​,
where β1β1​ captures the effect of ATP13A2 c.1306 LOF or CGP treatment, depending on the comparison. This full model was compared to a reduced null model lacking the treatment term (i.e., log⁡(πi/(1−πi))=β0log(πi​/(1−πi​))=β0​) to test whether inclusion of the treatment significantly improved model fit. Hypothesis tests were based on Wald or likelihood-ratio statistics, with p-values subsequently adjusted for multiple testing using a false discovery rate (FDR) procedure.
Unless otherwise specified, differentially methylated cytosines or regions (DMRs) were defined as those with FDR-adjusted q-values < 0.01 and a minimum absolute methylation difference exceeding a specified threshold. For stringent analyses focusing on the most robust changes, a methylation difference cutoff of >75 percentage points was applied. For some comparisons (e.g., exploratory analysis of c.1306 untreated vs c.1306+CGP), a more permissive cutoff of >25 percentage points was used to capture potentially treatment-responsive regions. In all cases, regions with higher methylation in the experimental group relative to the reference group were classified as hypermethylated, and those with lower methylation were classified as hypomethylated.
Comparative analyses across experimental conditions
Differential methylation analyses were carried out for multiple pairwise contrasts. The primary comparison examined DMRs between ATP13A2 c.1306 astrocytes and WT isogenic controls to identify methylation changes associated with ATP13A2 LOF. A secondary comparison evaluated c.1306 astrocytes treated with CGP versus untreated c.1306 cells to assess whether CGP exposure modulated methylation patterns within the mutant background. For each contrast, calculateDiffMeth() was applied to the relevant subset of samples, and the same statistical framework and FDR correction were used. The sign of the estimated methylation difference (experimental minus reference) was tracked for all regions to distinguish hypermethylated from hypomethylated DMRs.
Gene-level annotation of differentially methylated regions
To relate DMRs to genes and genomic features, significant regions were annotated using the UCSC Genome Browser. Genomic coordinates (chromosome, start, end, and strand) for DMRs were uploaded to the UCSC Table Browser and mapped to the human genome assembly hg38. Canonical gene annotations were obtained using the “knownCanonical” and related tables. For each DMR, the output included chromosome, start and end positions, gene symbol, and descriptive annotation fields (e.g. hg38 description).
Regions were assigned to genes based on overlap with gene bodies, promoters, or close proximity to annotated transcriptional units, yielding a gene-level list of loci associated with differential methylation. These annotation tables were used for integrating methylation data with gene expression results and for downstream enrichment analyses.
Integration of DNA methylation with gene expression data
To explore relationships between DNA methylation and transcriptional changes, DMR-associated genes were cross-referenced with differential gene expression results (e.g., log2 fold changes from RNA-seq or other expression profiling performed on matched or related samples). For each gene with associated hyper- or hypomethylated regions, methylation differences (percentage change between conditions) were paired with corresponding expression log2 fold-change values.
Scatterplots were generated to visualize the relationship between methylation difference and gene expression change in sets of hypermethylated and hypomethylated genes, and simple linear regression lines were fitted to summarize overall trends. In addition, heatmaps of the top hyper- and hypomethylated genes (ranked by methylation difference and/or significance) were constructed, displaying both methylation differences and expression log2 fold changes side-by-side. These integrative visualizations provided an exploratory view of how regional methylation changes might be associated with changes in gene expression, without imposing a priori assumptions about directionality.
Functional enrichment and pathway analysis
To identify biological functions and pathways associated with methylation changes, gene sets derived from DMR annotation were subjected to functional enrichment analysis. Lists of genes associated with hypermethylated and/or hypomethylated regions were uploaded to DAVID Functional Annotation Bioinformatics for over-representation analysis. DAVID was used to test for enrichment of gene ontology (GO) terms and curated functional categories relative to an appropriate background gene set.
Complementary enrichment analyses were performed in R using the clusterProfiler package. DMR-associated gene lists were input into clusterProfiler functions to test for over-representation of GO categories and other annotation sets supported by the package. Enrichment analyses were performed separately for different contrasts (e.g. c.1306 vs WT, c.1306+CGP vs c.1306) and for hyper- versus hypomethylated genes when relevant.
All enrichment analyses were conducted in an exploratory, hypothesis-generating framework. Only statistically supported terms (based on adjusted p- or q-value thresholds) were considered for interpretation, and results were integrated with knowledge of Parkinson’s disease biology and ATP13A2 function where appropriate.
Scope of the computational analysis
Collectively, these steps comprise a comprehensive computational workflow for RRBS data, including read preprocessing, bisulfite-aware alignment, CpG-level methylation calling, stringent QC, global and regional methylation quantification, statistically robust identification of differentially methylated regions with overdispersion correction, gene-level annotation, integration with gene expression data, and functional enrichment analyses. 

Specific details(If needed):


Experimental context
I analyzed RRBS-based DNA methylation data from human astrocytes to test whether ATP13A2 loss-of-function (c.1306) alters global and regional CpG methylation.
I framed this in the context of AMD1 upregulation, the SAM → dcSAM shift, and potential consequences for DNA/histone methylation relevant to Parkinson’s disease (PD).
I worked with three astrocyte conditions:
Wild-type (WT) isogenic control lines.
ATP13A2 c.1306 loss-of-function (LOF) mutant lines.
ATP13A2 c.1306 lines treated with CGP.
I analyzed approximately three biological replicates per condition.

RRBS preprocessing, alignment, and QC
I received raw RRBS FASTQ files for each astrocyte sample.
I trimmed each FASTQ file using CutAdapt:
I removed sequencing adapters.
I removed primer sequences.
I trimmed poly-A tails.
I removed any unwanted low-quality sequence from read ends.
After trimming, I performed post-trimming QC with FASTQC for every sample:
I checked per-base sequence quality and quality score distributions.
I examined GC content distributions.
I inspected sequence length distributions.
I checked for overrepresented sequences and residual adapter contamination.
I aligned trimmed reads to the human reference genome using Bismark:
I used Bismark’s bisulfite-aware alignment mode.
I extracted alignment efficiency and proportion of uniquely mapped reads.
I evaluated strand bias in the alignments.
I examined coverage across CpG sites from the Bismark reports.
I used Samtools on the Bismark BAM files:
I carried out alignment QC (basic mapping statistics, read counts, etc.).
I verified spike-in coverage to ensure appropriate representation of control sequences.
I used MethylDackel to assess bisulfite conversion and methylation calling:
I computed bisulfite conversion efficiency using non-CpG cytosines.
I measured bisulfite conversion efficiency using spike-in controls.
I generated methylation ratio estimates at CpG sites.


Methylation call files and import into R / methylKit
From the RRBS pipeline, I obtained CpG-level methylation call files (bedGraph/BAM-derived) via MultiQC.
For large methylation call files, I downloaded them sequentially using curl in Linux or Windows PowerShell to avoid transfer issues.
I imported the methylation call files into R using the methylKit package:
I used methylKit::methRead() with appropriate sample IDs, treatment labels, and genome assembly.
Each sample was read in as a MethylRaw object.
I combined all MethylRaw objects into a MethylRawList object for joint processing.
I confirmed that each CpG record included:
Chromosome.
Genomic position (start/end).
Strand.
Total coverage (read depth).
Counts of methylated vs unmethylated cytosines.

QC and descriptive methylation statistics
I summarized basic sequencing and methylation metrics for each sample:
Total read counts.
Number and fraction of uniquely aligned reads.
Bisulfite conversion rates from spike-ins and non-CpG sites.
Number of unique CpGs covered (~10–15 million CpGs per sample).
Average CpG depth (~11–15× per sample).
I used getMethylationStats() in methylKit:
To compute summary statistics for CpG methylation percentages.
To generate histograms of methylation distributions per sample.
I used getCoverageStats():
To compute summary statistics for CpG coverage.
To generate histograms and violin plots of coverage per sample.
I generated and inspected histograms and violin plots for:
The distribution of CpG methylation percentages.
The distribution of CpG coverage.
Based on these metrics and plots, I confirmed:
Bisulfite conversion efficiency was >99% for non-CpG sites and spike-ins.
Coverage and methylation distributions were appropriate and comparable across all samples.
All samples passed QC and were retained for downstream analyses.

Building the unified methylation object and global methylation analysis
I merged per-sample methylation data into a joint MethylBase object using unite() in methylKit:
I kept CpG sites that met the required minimum coverage and presence thresholds across samples.
The resulting MethylBase object contained, for each CpG:
Coverage per sample.
Methylated and unmethylated counts per sample.
I computed CpG percent methylation per sample with percMethylation().
I calculated global methylation levels by taking the mean CpG methylation per sample via colMeans().
I compared global CpG methylation between:
WT vs ATP13A2 c.1306 astrocytes.
(Optionally) c.1306 untreated vs c.1306+CGP.
I observed that ATP13A2 c.1306 astrocytes had ~0.4% lower global CpG methylation than WT and treated this as a descriptive trend, not a definitive epigenetic reprogramming result.

Regional methylation analysis using genomic tiling
To study regional methylation patterns, I applied tileMethylCounts() in methylKit:
I tiled the genome into fixed-width windows (genomic tiles).
I aggregated CpG methylated and unmethylated counts within each tile for each sample.
I converted the tiled methylation data into a MethylBase object representing genomic tiles rather than single CpGs.
For each tile, I had:
Tile coordinates (chr, start, end, strand).
Coverage (total counts).
Methylated vs unmethylated counts per sample.
Tile-level percent methylation per sample.

Differential methylation modeling and overdispersion correction
I used calculateDiffMeth() in methylKit to identify differentially methylated CpGs or tiles (DMRs).
For each site or tile, I modeled the number of methylated reads as:
methi∼Binomial(ni,πi)methi​∼Binomial(ni​,πi​),
where nini​ is the coverage and πiπi​ is the methylation proportion.
I treated the expected mean and variance as:
μi=niπiμi​=ni​πi​,
Var(methi)=niπi(1−πi)Var(methi​)=ni​πi​(1−πi​).
I accounted for overdispersion (variance > binomial expectation):
I estimated a scaling/dispersion parameter ϕ=χ2/(N−P)ϕ=χ2/(N−P), where:
NN = number of observations,
PP = number of model parameters.
I used this dispersion factor to inflate variance estimates and obtain more robust standard errors and p-values.
I fitted a logistic regression model per site/tile:
Full model:
log⁡(πi/(1−πi))=β0+β1Treatmentilog(πi​/(1−πi​))=β0​+β1​Treatmenti​.
Null model:
log⁡(πi/(1−πi))=β0log(πi​/(1−πi​))=β0​.
I tested whether including Treatment significantly improved model fit:
I obtained p-values for the treatment effect.
I corrected all p-values for multiple testing using FDR (q-values).
I defined differentially methylated sites/regions using:
q-value < 0.01, and
minimum absolute methylation difference thresholds:
Usually > 75% for very stringent DMRs.
For some exploratory treatment comparisons (e.g., c.1306 vs c.1306+CGP), I also considered regions with methylation difference > 25%.
I classified each significant region as:
Hypermethylated (higher methylation in experimental vs reference condition).
Hypomethylated (lower methylation in experimental vs reference condition).

Comparative DMR analyses across conditions
I ran separate calculateDiffMeth() analyses for:
ATP13A2 c.1306 vs WT astrocytes.
c.1306+CGP vs c.1306 untreated astrocytes.
In each comparison, I:
Used the same FDR threshold (q < 0.01).
Applied the relevant methylation difference cutoff (e.g., >75% or >25%).
Tracked the direction of methylation change (hyper vs hypo).
I summarized and visualized DMRs by:
Listing representative regions with chr/start/end, p-values, q-values, and methylation differences.
Plotting methylation difference vs –log10(q-value) to distinguish hyper- and hypomethylated regions.
Calculating the percentage of hyper vs hypo DMRs per chromosome.

Annotation of DMRs to genes
I exported significant DMR coordinates (chr, start, end, strand) to text files.
I uploaded these DMR coordinate files to the UCSC Genome Browser – Table Browser.
I used the human genome assembly hg38 and canonical gene annotation tables (e.g., knownCanonical):
I selected relevant fields such as chromosome, start, end, gene symbol, and description (hg38_desc).
I downloaded the annotation output and merged it with my DMR tables.
For each DMR, I assigned:
The corresponding gene symbol (where overlapping or nearest).
Additional annotation fields for downstream analysis.
I created gene-level lists of:
Genes associated with hypermethylated regions.
Genes associated with hypomethylated regions.
Genes showing methylation changes in the CGP treatment contrast.

Integration with gene expression and visualization
For genes associated with DMRs, I integrated methylation differences with differential expression (log2 fold change) values (from RNA-seq or other DEG data available for the same system).
I generated scatterplots of:
Gene expression (log2 fold change) vs methylation difference for hypermethylated genes.
Gene expression (log2 fold change) vs methylation difference for hypomethylated genes.
I overlaid linear regression trend lines on these scatterplots to visualize overall relationships between methylation change and expression change.
I created a heatmap of the top 10 hypermethylated and top 10 hypomethylated genes, combining:
Methylation differences.
Corresponding log2 fold changes in expression.
I used these plots as exploratory, hypothesis-generating visualizations rather than definitive mechanistic proof.

Functional enrichment and pathway analysis
I took the DMR-associated gene lists (e.g., hyper and hypo sets from each comparison) and performed pathway/GO enrichment analyses.
I used DAVID Functional Annotation Bioinformatics:
I uploaded gene lists with official gene symbols.
I ran over-representation analysis to identify enriched GO terms and functional categories.
I also used clusterProfiler in R:
I imported the same gene lists.
I ran GO and pathway enrichment (e.g., enrichGO, related functions).
I summarized enriched terms based on adjusted p- or q-value cutoffs.
I treated all enrichment analyses as exploratory:
I used the results to highlight pathways and processes potentially linked to ATP13A2 LOF–associated methylation changes and PD-relevant biology.


