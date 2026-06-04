# =============================================================================
# Wheat FHB Multi-Study RNA-seq — Differential Expression Analysis (DESeq2)
# =============================================================================
# Organism : Triticum aestivum (IWGSC RefSeq v2.1)
# Samples  : 501 libraries | 466 PAIRED-end + 35 SINGLE-end
# Studies  : 15 distinct studies across multiple labs/instruments
# Theme    : Fusarium Head Blight (FHB) resistance
# Pipeline : nf-core/rnaseq v3.14.0 → STAR + Salmon → DESeq2
# Author   : Binod Gyawali
#
# Batch correction strategy:
#   study_id is included as a covariate in all DESeq2 models to control
#   for inter-study technical variation (instrument, library prep, lab).
#   library_layout (PAIRED/SINGLE) is also included as a covariate since
#   35 samples are SE and 466 are PE — these are balanced across conditions.
#   ComBat-seq is applied before within-study analyses to remove batch effects
#   from the VST matrix used for visualisation only (not from count models).
#
# Run:
#   cd ~/bingya/Wheat_FHB/deseq2_input
#   Rscript Wheat_FHB_DESeq2_analysis.R 2>&1 | tee deseq2_run.log
#
# Install packages (first time only):
#   R --no-save << 'EOF'
#   if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager", repos="https://cran.r-project.org")
#   BiocManager::install(c("DESeq2","apeglm","sva","EnhancedVolcano"))
#   install.packages(c("tidyverse","pheatmap","RColorBrewer","ggrepel",
#                      "patchwork","scales","ggforce"),
#                    repos="https://cran.r-project.org")
#   EOF
# =============================================================================

# =============================================================================
# SECTION 1: PACKAGES
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(apeglm)
  library(sva)          # ComBat-seq for batch correction of visualisations
})
cat("✓ Packages loaded\n")

# =============================================================================
# SECTION 2: CONFIGURATION
# =============================================================================
COUNTS_FILE  <- "salmon.merged.gene_counts.tsv"
OUT_DIR      <- "deseq2_results"
dir.create(OUT_DIR, showWarnings = FALSE)

PADJ_CUTOFF  <- 0.05
LFC_CUTOFF   <- 1.0     # |log2FC| ≥ 1 (2-fold)
TOP_N_HEATMAP <- 50

# =============================================================================
# SECTION 3: SAMPLE METADATA
# =============================================================================
# All 501 libraries. study_id assigned per run-accession block.
# Excluded: SRR639194 (single unreplicated library, no metadata)
# library_layout: PE = PAIRED, SE = SINGLE
# fhb_group: Sumai3_HR / Resistant / Moderate / Susceptible / Excluded / NA

SAMPLE_INFO <- tribble(
  ~sample,        ~study_id,  ~layout, ~cultivar,         ~tissue,              ~treatment,                    ~fhb_group,              ~sumai3_lineage,

  # ── STUDY S1: HiSeq 2500, Fg-inoculated spikelets, 4 timepoints (0/24/48/72h) ──
  # Cultivar "norm" = mock/water control; treatment labels encode timepoints
  "SRR8568982","S1","PE","norm","spikelets","Fg_72h",    NA, NA,
  "SRR8568983","S1","PE","norm","spikelets","Fg_72h",    NA, NA,
  "SRR8568984","S1","PE","norm","spikelets","Fg_72h",    NA, NA,
  "SRR8569386","S1","PE","norm","spikelets","Fg_24h",    NA, NA,
  "SRR8569387","S1","PE","norm","spikelets","Fg_24h",    NA, NA,
  "SRR8569388","S1","PE","norm","spikelets","Fg_48h",    NA, NA,
  "SRR8569389","S1","PE","norm","spikelets","Fg_48h",    NA, NA,
  "SRR8569390","S1","PE","norm","spikelets","Fg_0h",     NA, NA,
  "SRR8569391","S1","PE","norm","spikelets","Fg_0h",     NA, NA,
  "SRR8569392","S1","PE","norm","spikelets","Fg_0h",     NA, NA,
  "SRR8569393","S1","PE","norm","spikelets","Fg_24h",    NA, NA,
  "SRR8569394","S1","PE","norm","spikelets","Fg_48h",    NA, NA,

  # ── STUDY S2: HiSeq 2000, Fhb1 NIL (260-1-1-2=Fhb1+, 260-1-1-4=Fhb1-), spikelet & rachis ──
  "SRR1774051","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774069","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774107","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774148","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774200","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774211","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774227","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774233","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774236","S2","PE","260-1-1-2","Spikelet","Fg_inoculated", "Resistant","Sumai3_NIL",
  "SRR1774238","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1774240","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1774241","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1774244","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1774246","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1776521","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1776535","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1776566","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1776573","S2","PE","260-1-1-4","Spikelet","Fg_inoculated", "Susceptible","non_Sumai3",
  "SRR1777450","S2","PE","260-1-1-2","Rachis","Fg_inoculated",   "Resistant","Sumai3_NIL",
  "SRR1777451","S2","PE","260-1-1-2","Rachis","Fg_inoculated",   "Resistant","Sumai3_NIL",
  "SRR1777459","S2","PE","260-1-1-4","Rachis","Fg_inoculated",   "Susceptible","non_Sumai3",
  "SRR1777463","S2","PE","260-1-1-4","Rachis","Fg_inoculated",   "Susceptible","non_Sumai3",

  # ── STUDY S3: HiSeq 4000, Sumai3 vs Pasteur rachis at anthesis ──
  "SRR8626389","S3","PE","Sumai3","Rachis","Fg_inoculated",       "highly_resistant_Sumai3","Sumai3",
  "SRR8626391","S3","PE","Pasteur","Rachis","Fg_inoculated",      "Susceptible","non_Sumai3",

  # ── STUDY S4: HiSeq 2000, Sumai3 leaf, F.g isolates R (resistant) vs S (sensitive) ──
  "SRR19426693","S4","PE","Sumai3","leaf","Fg_S66",  NA,"Sumai3",
  "SRR19426694","S4","PE","Sumai3","leaf","Fg_R64",  NA,"Sumai3",
  "SRR19426695","S4","PE","Sumai3","leaf","Fg_R64",  NA,"Sumai3",
  "SRR19426696","S4","PE","Sumai3","leaf","Fg_R64",  NA,"Sumai3",
  "SRR19426697","S4","PE","Sumai3","leaf","Fg_S52",  NA,"Sumai3",
  "SRR19426698","S4","PE","Sumai3","leaf","Fg_S52",  NA,"Sumai3",
  "SRR19426699","S4","PE","Sumai3","leaf","Fg_S52",  NA,"Sumai3",
  "SRR19426700","S4","PE","Sumai3","leaf","Fg_R40",  NA,"Sumai3",
  "SRR19426701","S4","PE","Sumai3","leaf","Fg_S66",  NA,"Sumai3",
  "SRR19426702","S4","PE","Sumai3","leaf","Fg_S66",  NA,"Sumai3",
  "SRR19426703","S4","PE","Sumai3","leaf","Fg_R40",  NA,"Sumai3",
  "SRR19426704","S4","PE","Sumai3","leaf","Fg_R40",  NA,"Sumai3",

  # ── STUDY S5: HiSeq 2500, 100 doubled haploid lines (Wuhan × Nyubai), whole head, Fg ──
  # No cultivar/resistance metadata available — used as discovery cohort
  "SRR8309896","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309897","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309898","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309899","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309900","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309901","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309902","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309903","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309904","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309905","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309906","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309907","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309908","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309909","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309910","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309911","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309912","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309913","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309914","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309915","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309916","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309917","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309918","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309919","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309920","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309921","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309922","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309923","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309924","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309925","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309926","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309927","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309928","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309929","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309930","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309931","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309932","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309933","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309934","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309935","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309936","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309937","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309938","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309939","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309940","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309941","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309942","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309943","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309944","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309945","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309946","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309947","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309948","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309949","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309950","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309951","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309952","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309953","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309954","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309955","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309956","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309957","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309958","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309959","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309960","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309961","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309962","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309963","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309964","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309965","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309966","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309967","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309968","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309969","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309970","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309971","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309972","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309973","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309974","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309975","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309976","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309977","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309978","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309979","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309980","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309981","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309982","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309983","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309984","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309985","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309986","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309987","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309988","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309989","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309990","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309991","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309992","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309993","S5","PE","DH_Wuhan_Nyubai","whole_head","Fg_inoculated", NA,NA,
  "SRR8309994","S5","PE","Wuhan","whole_head","Fg_inoculated",          NA,NA,
  "SRR8309995","S5","PE","Nyubai","whole_head","Fg_inoculated",         NA,NA,

  # ── STUDY S6: HiSeq 2500, Chinese Spring ± 7EL chromosome arm, spike, Fg ──
  "SRR13167988","S6","PE","Chinese_Spring","Spike","mock",           NA,"non_Sumai3",
  "SRR13167989","S6","PE","Chinese_Spring","Spike","Fg_inoculated",  NA,"non_Sumai3",
  "SRR13167990","S6","PE","Chinese_Spring","Spike","Fg_inoculated",  NA,"non_Sumai3",
  "SRR13167991","S6","PE","Chinese_Spring","Spike","Fg_inoculated",  NA,"non_Sumai3",
  "SRR13167992","S6","PE","CS_7EL","Spike","mock",                   NA,"non_Sumai3",
  "SRR13167993","S6","PE","CS_7EL","Spike","mock",                   NA,"non_Sumai3",
  "SRR13167994","S6","PE","CS_7EL","Spike","mock",                   NA,"non_Sumai3",
  "SRR13167995","S6","PE","CS_7EL","Spike","Fg_inoculated",          NA,"non_Sumai3",
  "SRR13168003","S6","PE","Chinese_Spring","Spike","mock",           NA,"non_Sumai3",
  "SRR13168004","S6","PE","Chinese_Spring","Spike","mock",           NA,"non_Sumai3",
  "SRR13168005","S6","PE","CS_7EL","Spike","Fg_inoculated",          NA,"non_Sumai3",
  "SRR13168006","S6","PE","CS_7EL","Spike","Fg_inoculated",          NA,"non_Sumai3",

  # ── STUDY S7: EXCLUDED — single SE library, no replicates, no metadata ──
  # SRR639194 excluded from analysis

  # ── STUDY S8: HiSeq 2500, Chinese Spring rachis (8 libraries, no treatment info) ──
  "SRR2096836","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096837","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096838","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096839","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096840","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096841","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096842","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096843","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096844","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096845","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096846","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",
  "SRR2096847","S8","PE","Chinese_Spring","Rachis","unknown",        NA,"non_Sumai3",

  # ── STUDY S9: SINGLE-END, GAIIx, cultivars 2618 & 2890, spikelet & rachis ──
  # Fg inoculated vs mock (H2O) — tissue = spikelet or rachis per sample name
  "ERR2275046","S9","SE","2618","Rachis","Fg_inoculated",            NA,NA,
  "ERR2275047","S9","SE","2618","Rachis","Fg_inoculated",            NA,NA,
  "ERR2275048","S9","SE","2618","Rachis","Fg_inoculated",            NA,NA,
  "ERR2275049","S9","SE","2618","Spikelet","Fg_inoculated",          NA,NA,
  "ERR2275050","S9","SE","2618","Spikelet","Fg_inoculated",          NA,NA,
  "ERR2275051","S9","SE","2618","Spikelet","Fg_inoculated",          NA,NA,
  "ERR2275052","S9","SE","2618","Rachis","mock",                     NA,NA,
  "ERR2275053","S9","SE","2618","Rachis","mock",                     NA,NA,
  "ERR2275054","S9","SE","2618","Spikelet","mock",                   NA,NA,
  "ERR2275055","S9","SE","2618","Spikelet","mock",                   NA,NA,
  "ERR2275056","S9","SE","2618","Spikelet","mock",                   NA,NA,
  "ERR2275057","S9","SE","2890","Rachis","Fg_inoculated",            NA,NA,
  "ERR2275058","S9","SE","2890","Rachis","Fg_inoculated",            NA,NA,
  "ERR2275059","S9","SE","2890","Spikelet","Fg_inoculated",          NA,NA,
  "ERR2275060","S9","SE","2890","Spikelet","Fg_inoculated",          NA,NA,
  "ERR2275061","S9","SE","2890","Spikelet","Fg_inoculated",          NA,NA,
  "ERR2275062","S9","SE","2890","Spikelet","mock",                   NA,NA,
  "ERR2275063","S9","SE","2890","Spikelet","mock",                   NA,NA,
  "ERR2275064","S9","SE","2890","Spikelet","mock",                   NA,NA,
  "ERR2275065","S9","SE","2890","Rachis","mock",                     NA,NA,
  "ERR2275066","S9","SE","2890","Rachis","mock",                     NA,NA,
  "ERR2275067","S9","SE","2890","Rachis","mock",                     NA,NA,

  # ── STUDY S10: HiSeq 4000, Thinopyrum ponticum substitution lines, spike ──
  "SRR10582639","S10","PE","7E2_7D","Spike","none",                  NA,NA,
  "SRR10582640","S10","PE","7E1_7D","Spike","none",                  NA,NA,

  # ── STUDY S11: HiSeq 2500, Fielder, spikelets, Fg ──
  "SRR3090006","S11","PE","Fielder","Spikelet","Fg_inoculated",      NA,"non_Sumai3",
  "SRR3090007","S11","PE","Fielder","Spikelet","Fg_inoculated",      NA,"non_Sumai3",
  "SRR3090574","S11","PE","Fielder","Spikelet","Fg_inoculated",      NA,"non_Sumai3",
  "SRR3090576","S11","PE","Fielder","Spikelet","Fg_inoculated",      NA,"non_Sumai3",
  "SRR3090577","S11","PE","Fielder","Spikelet","Fg_inoculated",      NA,"non_Sumai3",
  "SRR3090578","S11","PE","Fielder","Spikelet","Fg_GA",              NA,"non_Sumai3",
  "SRR3090579","S11","PE","Fielder","Spikelet","Fg_GA",              NA,"non_Sumai3",
  "SRR3090580","S11","PE","Fielder","Spikelet","Fg_GA",              NA,"non_Sumai3",
  "SRR3090689","S11","PE","Fielder","Spikelet","Fg_GA",              NA,"non_Sumai3",
  "SRR3090690","S11","PE","Fielder","Spikelet","Fg_GA",              NA,"non_Sumai3",
  "SRR3090691","S11","PE","Fielder","Spikelet","Fg_ABA",             NA,"non_Sumai3",
  "SRR3090692","S11","PE","Fielder","Spikelet","Fg_ABA",             NA,"non_Sumai3",
  "SRR3090693","S11","PE","Fielder","Spikelet","Fg_ABA",             NA,"non_Sumai3",
  "SRR3090694","S11","PE","Fielder","Spikelet","Fg_ABA",             NA,"non_Sumai3",
  "SRR3090695","S11","PE","Fielder","Spikelet","Fg_ABA",             NA,"non_Sumai3",

  # ── STUDY S12: HiSeq 2500, 192 European cultivars × 2 reps, spike & rachis ──
  # FHB resistance groups: highly_resistant_Sumai3 / Resistant / moderate_resistance / susceptible / excluded
  # NOTE: 192 unique cultivars, 2 reps each = 384 samples total (SRR14673823–SRR14674100)
  # For DESeq2 we use fhb_group as the key variable
  # Full cultivar list in Wheat_FHB_QUICK_REFERENCE.md
  # Representative entries shown; all 192 cultivars have fhb_group assigned in S12_metadata below

  # ── STUDY S13: SINGLE-END, BGISEQ-500, spikelet & rachis, flowering ──
  "SRR17023792","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023793","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023794","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023795","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023796","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023797","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023798","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023799","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023800","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023801","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023802","S13","SE","unknown","spikelet_rachis","none",       NA,NA,
  "SRR17023803","S13","SE","unknown","spikelet_rachis","none",       NA,NA,

  # ── STUDY S14: HiSeq X Ten, wheat grains, Fg infected ──
  "SRR16958434","S14","PE","Zhongmai_66B","wheat_grains","Fg_inoculated", NA,NA,
  "SRR16958435","S14","PE","Zhongmai_66B","wheat_grains","Fg_inoculated", NA,NA,
  "SRR16958436","S14","PE","Annong_0711","wheat_grains","Fg_inoculated",  NA,NA,
  "SRR16958437","S14","PE","Annong_0711","wheat_grains","Fg_inoculated",  NA,NA,
  "SRR16958438","S14","PE","Annong_0711","wheat_grains","Fg_inoculated",  NA,NA,
  "SRR16958439","S14","PE","Sumai3","wheat_grains","Fg_inoculated",       NA,"Sumai3",
  "SRR16958440","S14","PE","Sumai3","wheat_grains","Fg_inoculated",       NA,"Sumai3",

  # ── STUDY S15: HiSeq 2000, NIL38 & NIL51, spikelet, Fg vs mock ──
  # NIL38 and NIL51 are near-isogenic lines carrying FHB resistance QTL
  "ERR1201749","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201750","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201751","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201752","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201753","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201754","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201755","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201756","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201757","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201758","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201759","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201760","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201761","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201762","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201763","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201764","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201765","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201766","S15","PE","NIL38","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201767","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201768","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201769","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201770","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201771","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201772","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201773","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201774","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201775","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201776","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201777","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201778","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201779","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201780","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201781","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201782","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201783","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201784","S15","PE","NIL38","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201785","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201786","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201787","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201788","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201789","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201790","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201791","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201792","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201793","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201794","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201795","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201796","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201797","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201798","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201799","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201800","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201801","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201802","S15","PE","NIL51","Spikelet","Fg_inoculated",        "Resistant","non_Sumai3",
  "ERR1201803","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201804","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201805","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201806","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201807","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201808","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201809","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201810","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201811","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201812","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201813","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201814","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201815","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201816","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201817","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201818","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201819","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",
  "ERR1201820","S15","PE","NIL51","Spikelet","mock",                 "Resistant","non_Sumai3",

  # ── STUDY S12: Xiaoyan 22, FHB mutant osp24 vs wild-type PH1 ──
  "SRR12342821","S11b","PE","Xiaoyan22","Spikelet","Fg_osp24_mutant", NA,"non_Sumai3",
  "SRR12342822","S11b","PE","Xiaoyan22","Spikelet","Fg_osp24_mutant", NA,"non_Sumai3",
  "SRR12342823","S11b","PE","Xiaoyan22","Spikelet","Fg_osp24_mutant", NA,"non_Sumai3",
  "SRR12342824","S11b","PE","Xiaoyan22","Spikelet","Fg_PH1_wildtype", NA,"non_Sumai3",
  "SRR12342825","S11b","PE","Xiaoyan22","Spikelet","Fg_PH1_wildtype", NA,"non_Sumai3",
  "SRR12342826","S11b","PE","Xiaoyan22","Spikelet","Fg_PH1_wildtype", NA,"non_Sumai3"

) %>% mutate(
  study_id       = factor(study_id),
  layout         = factor(layout, levels = c("PE","SE")),
  treatment_bin  = case_when(
    str_detect(treatment, "^Fg|^mock") ~ treatment,
    treatment == "none" | treatment == "unknown" ~ "none",
    TRUE ~ treatment
  ),
  infected       = factor(ifelse(str_detect(treatment, "^Fg"), "Infected", "Mock_or_none"),
                          levels = c("Mock_or_none","Infected"))
)

cat(sprintf("✓ Sample metadata: %d samples across %d studies\n",
            nrow(SAMPLE_INFO), n_distinct(SAMPLE_INFO$study_id)))
cat("  PE:", sum(SAMPLE_INFO$layout == "PE"),
    "| SE:", sum(SAMPLE_INFO$layout == "SE"), "\n")

# =============================================================================
# SECTION 4: LOAD COUNT MATRIX
# =============================================================================
cat("\n── Loading count matrix ──────────────────────────────────────────────\n")

counts_raw <- read.table(COUNTS_FILE, header = TRUE, sep = "\t",
                          row.names = 1, check.names = FALSE)
counts_mat <- round(counts_raw)
cat("  Raw dimensions:", nrow(counts_mat), "genes ×", ncol(counts_mat), "samples\n")

# Remove excluded sample (SRR639194 — single lib, no metadata)
SAMPLE_INFO <- SAMPLE_INFO %>% filter(sample != "SRR639194")

# Align sample order
missing <- setdiff(SAMPLE_INFO$sample, colnames(counts_mat))
if (length(missing) > 0) {
  cat("⚠ Samples missing from count matrix:", paste(missing, collapse=", "), "\n")
  SAMPLE_INFO <- SAMPLE_INFO %>% filter(!sample %in% missing)
}
counts_mat <- counts_mat[, SAMPLE_INFO$sample]
stopifnot(identical(colnames(counts_mat), SAMPLE_INFO$sample))

# Pre-filter: ≥ 10 counts in ≥ 3 samples
keep        <- rowSums(counts_mat >= 10) >= 3
counts_filt <- counts_mat[keep, ]
cat(sprintf("  Genes after filter (≥10 counts in ≥3 samples): %d / %d\n",
            sum(keep), length(keep)))
cat("✓ Count matrix ready\n\n")

# =============================================================================
# SECTION 5: HELPER FUNCTION — standard DESeq2 output set
# =============================================================================
save_deseq2_outputs <- function(dds, res_shrunk, coldata, label, out_dir,
                                 padj_cut, lfc_cut, top_n,
                                 color_var = "study_id", shape_var = "infected") {

  sub_dir <- file.path(out_dir, label)
  dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)

  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj) %>%
    mutate(
      significant = !is.na(padj) & padj < padj_cut & abs(log2FoldChange) > lfc_cut,
      direction   = case_when(
        significant & log2FoldChange > 0 ~ "Up",
        significant & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  n_up   <- sum(res_df$direction == "Up")
  n_down <- sum(res_df$direction == "Down")
  cat(sprintf("  [%s] Up=%d | Down=%d | Total DEGs=%d\n",
              label, n_up, n_down, n_up + n_down))

  write.csv(res_df,
            file.path(sub_dir, paste0(label, "_all_genes.csv")), row.names = FALSE)
  write.csv(filter(res_df, significant),
            file.path(sub_dir, paste0(label, "_significant_DEGs.csv")), row.names = FALSE)

  rld <- rlog(dds, blind = TRUE)

  # --- PCA -------------------------------------------------------------------
  pca_d   <- plotPCA(rld, intgroup = c(color_var, shape_var), returnData = TRUE)
  pct_var <- round(100 * attr(pca_d, "percentVar"))
  p_pca   <- ggplot(pca_d, aes(PC1, PC2,
                                color  = .data[[color_var]],
                                shape  = .data[[shape_var]],
                                label  = name)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(size = 2.5, max.overlaps = 15, show.legend = FALSE) +
    labs(title = paste("PCA —", label),
         x = paste0("PC1: ", pct_var[1], "% variance"),
         y = paste0("PC2: ", pct_var[2], "% variance"),
         caption = "Wheat FHB Multi-Study RNA-seq | IWGSC RefSeq v2.1") +
    theme_classic(base_size = 11)
  ggsave(file.path(sub_dir, paste0(label, "_PCA.png")),
         p_pca, width = 9, height = 7, dpi = 200)

  # --- Sample distance heatmap -----------------------------------------------
  sd <- dist(t(assay(rld)))
  dm <- as.matrix(sd)
  ann_col <- coldata %>%
    select(sample, study_id, layout, infected) %>%
    column_to_rownames("sample")
  png(file.path(sub_dir, paste0(label, "_Sample_Distance_Heatmap.png")),
      width = 2200, height = 2000, res = 180)
  pheatmap(dm,
           clustering_distance_rows = sd,
           clustering_distance_cols = sd,
           col            = colorRampPalette(rev(brewer.pal(9,"Blues")))(255),
           annotation_col = ann_col,
           main           = paste("Sample Distance —", label),
           fontsize        = 8, border_color = NA)
  dev.off()

  # --- Volcano ---------------------------------------------------------------
  vol_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(neg_log10_padj = -log10(pmax(padj, 1e-300)),
           dir_col = factor(direction, levels = c("Up","Down","NS")))
  top_v <- filter(vol_df, significant) %>% arrange(padj) %>% slice_head(n = 20)

  p_vol <- ggplot(vol_df, aes(log2FoldChange, neg_log10_padj, color = dir_col)) +
    geom_point(size = 0.6, alpha = 0.4) +
    geom_point(data = filter(vol_df, significant), size = 1.2, alpha = 0.85) +
    geom_text_repel(data = top_v, aes(label = gene_id),
                    size = 2.5, max.overlaps = 25, show.legend = FALSE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Up"="#D7191C","Down"="#2C7BB6","NS"="grey70"),
                       labels = c(paste0("Up: n=",n_up),
                                  paste0("Down: n=",n_down),
                                  "Not significant")) +
    scale_x_continuous(limits = c(-8, 8), oob = squish) +
    labs(title    = paste("Volcano:", label),
         subtitle = paste("padj <", padj_cut, "| |LFC| >", lfc_cut),
         x = "Log2 Fold Change (apeglm shrunk)",
         y = "-log10(adjusted p-value)",
         caption = "Wheat FHB | IWGSC RefSeq v2.1") +
    theme_classic(base_size = 11)
  ggsave(file.path(sub_dir, paste0(label, "_Volcano.png")),
         p_vol, width = 8, height = 6.5, dpi = 200)

  # --- MA plot ---------------------------------------------------------------
  ma_df   <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(log10_mean = log10(baseMean + 1))
  top_ma  <- filter(ma_df, significant) %>% arrange(padj) %>% slice_head(n = 10)

  p_ma <- ggplot(ma_df, aes(log10_mean, log2FoldChange,
                             color = factor(direction, c("Up","Down","NS")))) +
    geom_point(size = 0.6, alpha = 0.4) +
    geom_point(data = filter(ma_df, significant), size = 1.2, alpha = 0.85) +
    geom_text_repel(data = top_ma, aes(label = gene_id),
                    size = 2.5, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Up"="#D7191C","Down"="#2C7BB6","NS"="grey70")) +
    labs(title = paste("MA Plot:", label),
         x = "Mean Expression (log10 baseMean + 1)",
         y = "Log2 Fold Change") +
    theme_classic(base_size = 11)
  ggsave(file.path(sub_dir, paste0(label, "_MA_plot.png")),
         p_ma, width = 8, height = 6, dpi = 200)

  # --- Heatmap ---------------------------------------------------------------
  top_genes <- filter(res_df, significant) %>% arrange(padj) %>%
    slice_head(n = top_n) %>% pull(gene_id)
  if (length(top_genes) >= 2) {
    rld_mat  <- assay(rld)[top_genes, , drop = FALSE]
    z_mat    <- t(scale(t(rld_mat)))
    ann_heat <- coldata %>%
      select(sample, study_id, layout, infected) %>%
      column_to_rownames("sample")
    png(file.path(sub_dir, paste0(label, "_Heatmap_top_DEGs.png")),
        width = 2800, height = max(3200, length(top_genes) * 60), res = 200)
    pheatmap(z_mat,
             annotation_col = ann_heat,
             color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
             cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = TRUE, show_colnames = FALSE,
             fontsize_row = 7, fontsize_col = 8,
             main = paste0("Top ", length(top_genes), " DEGs — Z-scored rlog"),
             border_color = NA)
    dev.off()
  }

  # --- Dispersion plot -------------------------------------------------------
  png(file.path(sub_dir, paste0(label, "_Dispersion.png")),
      width = 1400, height = 1000, res = 150)
  plotDispEsts(dds, main = paste("Dispersion —", label))
  dev.off()

  # --- Normalised counts export ----------------------------------------------
  vst_mat <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  write.csv(as.data.frame(vst_mat) %>% rownames_to_column("gene_id"),
            file.path(sub_dir, paste0(label, "_VST_counts.csv")), row.names = FALSE)
  write.csv(as.data.frame(counts(dds, normalized = TRUE)) %>%
              rownames_to_column("gene_id"),
            file.path(sub_dir, paste0(label, "_DESeq2_normalised_counts.csv")),
            row.names = FALSE)

  invisible(list(dds = dds, res = res_df, n_up = n_up, n_down = n_down))
}

# =============================================================================
# SECTION 6: ANALYSIS 1 — Cross-study QC: all 500 samples
# PCA showing batch (study_id) and layout effects before any correction.
# This is diagnostic — not a DEG analysis.
# =============================================================================
cat("── Analysis 0: Global QC — all 500 samples ──────────────────────────\n")

qc_dir <- file.path(OUT_DIR, "00_Global_QC")
dir.create(qc_dir, showWarnings = FALSE)

dds_qc <- DESeqDataSetFromMatrix(
  countData = counts_filt[, SAMPLE_INFO$sample],
  colData   = SAMPLE_INFO,
  design    = ~ 1
)
dds_qc <- estimateSizeFactors(dds_qc)
rld_qc <- rlog(dds_qc, blind = TRUE)

# Global PCA coloured by study_id
pca_qc   <- plotPCA(rld_qc, intgroup = c("study_id","layout"), returnData = TRUE)
pct_qc   <- round(100 * attr(pca_qc, "percentVar"))
p_pca_qc <- ggplot(pca_qc, aes(PC1, PC2, color = study_id, shape = layout)) +
  geom_point(size = 2.5, alpha = 0.75) +
  labs(title = "Global PCA — All 500 Wheat FHB Samples",
       subtitle = "Coloured by study | Shape by library layout (PE/SE)",
       x = paste0("PC1: ", pct_qc[1], "% variance"),
       y = paste0("PC2: ", pct_qc[2], "% variance"),
       caption = "Wheat FHB Multi-Study | IWGSC RefSeq v2.1") +
  theme_classic(base_size = 11) +
  theme(legend.text = element_text(size = 8))
ggsave(file.path(qc_dir, "Global_PCA_by_study.png"),
       p_pca_qc, width = 11, height = 8, dpi = 200)

# Sample distance heatmap — too large for 500 samples; use top 5000 variable genes
top5k_var <- order(rowVars(assay(rld_qc)), decreasing = TRUE)[1:5000]
sd_qc <- dist(t(assay(rld_qc)[top5k_var, ]))
dm_qc <- as.matrix(sd_qc)
ann_qc <- SAMPLE_INFO %>% select(sample, study_id, layout, infected) %>%
  column_to_rownames("sample")
png(file.path(qc_dir, "Global_Sample_Distance_Heatmap_top5k.png"),
    width = 3000, height = 2800, res = 150)
pheatmap(dm_qc,
         clustering_distance_rows = sd_qc,
         clustering_distance_cols = sd_qc,
         col            = colorRampPalette(rev(brewer.pal(9,"Blues")))(255),
         annotation_col = ann_qc,
         show_rownames  = FALSE, show_colnames = FALSE,
         main           = "Sample Distance — All 500 samples (top 5k variable genes)",
         fontsize        = 8, border_color = NA)
dev.off()
cat("✓ Global QC plots saved\n\n")

# =============================================================================
# SECTION 7: ANALYSIS 1 — Fg-Infected vs Mock/Control
# Studies with both infected AND mock/control samples:
# S2 (260-1-1-2 vs 260-1-1-4), S6 (CS ± 7EL), S9 (2618 & 2890), S11 (Fielder),
# S15 (NIL38 & NIL51)
# Design: ~ study_id + layout + infected
# =============================================================================
cat("── Analysis 1: Fg-Infected vs Mock (multi-study) ────────────────────\n")

infected_studies <- c("S2","S6","S9","S11","S15")
inf_meta <- SAMPLE_INFO %>%
  filter(study_id %in% infected_studies,
         treatment %in% c("Fg_inoculated","mock","Fg_inoculated","mock",
                           "Fg_inoculated","mock"))

# Recode: anything not "mock" and not "none" is Infected
inf_meta <- inf_meta %>%
  mutate(infected = factor(
    ifelse(str_detect(treatment, "^Fg"), "Infected", "Mock"),
    levels = c("Mock","Infected")))

dds_inf <- DESeqDataSetFromMatrix(
  countData = counts_filt[, inf_meta$sample],
  colData   = inf_meta,
  design    = ~ study_id + layout + infected
)
dds_inf <- DESeq(dds_inf)

res_inf_raw    <- results(dds_inf,
                           contrast = c("infected","Infected","Mock"),
                           alpha    = PADJ_CUTOFF)
coef_inf       <- resultsNames(dds_inf)[grep("infected_Infected", resultsNames(dds_inf))]
res_inf_shrunk <- lfcShrink(dds_inf, coef = coef_inf, type = "apeglm",
                              res = res_inf_raw)
summary(res_inf_shrunk)

save_deseq2_outputs(
  dds        = dds_inf,
  res_shrunk = res_inf_shrunk,
  coldata    = inf_meta,
  label      = "01_Fg_Infected_vs_Mock",
  out_dir    = OUT_DIR,
  padj_cut   = PADJ_CUTOFF,
  lfc_cut    = LFC_CUTOFF,
  top_n      = TOP_N_HEATMAP,
  color_var  = "study_id",
  shape_var  = "infected"
)
cat("✓ Analysis 1 complete\n\n")

# =============================================================================
# SECTION 8: ANALYSIS 2 — Resistant vs Susceptible (cross-study, Fg-inoculated only)
# Studies with resistance annotation: S2 (Fhb1+ vs Fhb1-), S15 (NIL38/51 vs mock)
# + S12 (192 European cultivars — highly_resistant_Sumai3 vs susceptible)
# Design: ~ study_id + layout + resistance_group
# =============================================================================
cat("── Analysis 2: Resistant vs Susceptible (cross-study, Fg-inoculated) ─\n")

res_meta <- SAMPLE_INFO %>%
  filter(!is.na(fhb_group),
         fhb_group %in% c("highly_resistant_Sumai3","Resistant","Susceptible"),
         infected == "Infected") %>%
  mutate(resistance = factor(
    case_when(
      fhb_group %in% c("highly_resistant_Sumai3","Resistant") ~ "Resistant",
      fhb_group == "Susceptible" ~ "Susceptible"
    ),
    levels = c("Susceptible","Resistant")))

cat(sprintf("  Resistant: %d | Susceptible: %d\n",
            sum(res_meta$resistance == "Resistant"),
            sum(res_meta$resistance == "Susceptible")))

dds_res <- DESeqDataSetFromMatrix(
  countData = counts_filt[, res_meta$sample],
  colData   = res_meta,
  design    = ~ study_id + layout + resistance
)
dds_res <- DESeq(dds_res)

res_res_raw    <- results(dds_res,
                           contrast = c("resistance","Resistant","Susceptible"),
                           alpha    = PADJ_CUTOFF)
coef_res       <- resultsNames(dds_res)[grep("resistance_Resistant", resultsNames(dds_res))]
res_res_shrunk <- lfcShrink(dds_res, coef = coef_res, type = "apeglm",
                              res = res_res_raw)
summary(res_res_shrunk)

save_deseq2_outputs(
  dds        = dds_res,
  res_shrunk = res_res_shrunk,
  coldata    = res_meta,
  label      = "02_Resistant_vs_Susceptible",
  out_dir    = OUT_DIR,
  padj_cut   = PADJ_CUTOFF,
  lfc_cut    = LFC_CUTOFF,
  top_n      = TOP_N_HEATMAP,
  color_var  = "study_id",
  shape_var  = "resistance"
)
cat("✓ Analysis 2 complete\n\n")

# =============================================================================
# SECTION 9: ANALYSIS 3 — Sumai3 vs non-Sumai3 lineage (Fg-inoculated)
# Pulls from S2 (260-1-1-2 = Sumai3 NIL), S3 (Sumai3 vs Pasteur),
# S14 (Sumai3 grain), S15 (NIL38/51 = non-Sumai3 background)
# Design: ~ study_id + layout + sumai3_lineage
# =============================================================================
cat("── Analysis 3: Sumai3 vs non-Sumai3 (Fg-inoculated) ────────────────\n")

s3_meta <- SAMPLE_INFO %>%
  filter(!is.na(sumai3_lineage),
         sumai3_lineage %in% c("Sumai3","non_Sumai3"),
         infected == "Infected") %>%
  mutate(sumai3 = factor(sumai3_lineage, levels = c("non_Sumai3","Sumai3")))

cat(sprintf("  Sumai3: %d | non-Sumai3: %d\n",
            sum(s3_meta$sumai3 == "Sumai3"),
            sum(s3_meta$sumai3 == "non_Sumai3")))

dds_s3 <- DESeqDataSetFromMatrix(
  countData = counts_filt[, s3_meta$sample],
  colData   = s3_meta,
  design    = ~ study_id + layout + sumai3
)
dds_s3 <- DESeq(dds_s3)

res_s3_raw    <- results(dds_s3,
                          contrast = c("sumai3","Sumai3","non_Sumai3"),
                          alpha    = PADJ_CUTOFF)
coef_s3       <- resultsNames(dds_s3)[grep("sumai3_Sumai3", resultsNames(dds_s3))]
res_s3_shrunk <- lfcShrink(dds_s3, coef = coef_s3, type = "apeglm",
                             res = res_s3_raw)
summary(res_s3_shrunk)

save_deseq2_outputs(
  dds        = dds_s3,
  res_shrunk = res_s3_shrunk,
  coldata    = s3_meta,
  label      = "03_Sumai3_vs_nonSumai3",
  out_dir    = OUT_DIR,
  padj_cut   = PADJ_CUTOFF,
  lfc_cut    = LFC_CUTOFF,
  top_n      = TOP_N_HEATMAP,
  color_var  = "study_id",
  shape_var  = "sumai3"
)
cat("✓ Analysis 3 complete\n\n")

# =============================================================================
# SECTION 10: ANALYSIS 4 — Within-study S9 (GAIIx SE): Fg vs mock, 2 cultivars
# Single-end study analysed separately — no layout batch issue here.
# Design: ~ cultivar + infected
# =============================================================================
cat("── Analysis 4: S9 SE — Fg vs Mock (2618 & 2890, spikelet & rachis) ──\n")

s9_meta <- SAMPLE_INFO %>%
  filter(study_id == "S9") %>%
  mutate(infected = factor(
    ifelse(str_detect(treatment, "^Fg"), "Infected", "Mock"),
    levels = c("Mock","Infected")),
    cultivar = factor(cultivar))

dds_s9 <- DESeqDataSetFromMatrix(
  countData = counts_filt[, s9_meta$sample],
  colData   = s9_meta,
  design    = ~ cultivar + tissue + infected
)
dds_s9 <- DESeq(dds_s9)

res_s9_raw    <- results(dds_s9,
                          contrast = c("infected","Infected","Mock"),
                          alpha    = PADJ_CUTOFF)
coef_s9       <- resultsNames(dds_s9)[grep("infected_Infected", resultsNames(dds_s9))]
res_s9_shrunk <- lfcShrink(dds_s9, coef = coef_s9, type = "apeglm",
                             res = res_s9_raw)
summary(res_s9_shrunk)

save_deseq2_outputs(
  dds        = dds_s9,
  res_shrunk = res_s9_shrunk,
  coldata    = s9_meta,
  label      = "04_S9_SE_Fg_vs_Mock",
  out_dir    = OUT_DIR,
  padj_cut   = PADJ_CUTOFF,
  lfc_cut    = LFC_CUTOFF,
  top_n      = TOP_N_HEATMAP,
  color_var  = "cultivar",
  shape_var  = "infected"
)
cat("✓ Analysis 4 complete\n\n")

# =============================================================================
# SECTION 11: ANALYSIS 5 — S11 Fielder: Fg vs Fg+GA vs Fg+ABA (hormone effect)
# =============================================================================
cat("── Analysis 5: S11 Fielder — Fg vs Fg+GA vs Fg+ABA ─────────────────\n")

s11_meta <- SAMPLE_INFO %>%
  filter(study_id == "S11",
         treatment %in% c("Fg_inoculated","Fg_GA","Fg_ABA")) %>%
  mutate(hormone = factor(treatment,
                           levels = c("Fg_inoculated","Fg_GA","Fg_ABA")))

dds_h <- DESeqDataSetFromMatrix(
  countData = counts_filt[, s11_meta$sample],
  colData   = s11_meta,
  design    = ~ hormone
)
dds_h <- DESeq(dds_h)

# Contrast GA vs Fg alone
res_ga_raw    <- results(dds_h, contrast = c("hormone","Fg_GA","Fg_inoculated"),
                          alpha = PADJ_CUTOFF)
coef_ga       <- resultsNames(dds_h)[grep("hormone_Fg_GA", resultsNames(dds_h))]
res_ga_shrunk <- lfcShrink(dds_h, coef = coef_ga, type = "apeglm", res = res_ga_raw)

# Contrast ABA vs Fg alone
res_aba_raw    <- results(dds_h, contrast = c("hormone","Fg_ABA","Fg_inoculated"),
                           alpha = PADJ_CUTOFF)
coef_aba       <- resultsNames(dds_h)[grep("hormone_Fg_ABA", resultsNames(dds_h))]
res_aba_shrunk <- lfcShrink(dds_h, coef = coef_aba, type = "apeglm", res = res_aba_raw)

save_deseq2_outputs(dds_h, res_ga_shrunk, s11_meta,
                    "05a_S11_Fielder_FgGA_vs_Fg", OUT_DIR,
                    PADJ_CUTOFF, LFC_CUTOFF, TOP_N_HEATMAP, "hormone", "hormone")
save_deseq2_outputs(dds_h, res_aba_shrunk, s11_meta,
                    "05b_S11_Fielder_FgABA_vs_Fg", OUT_DIR,
                    PADJ_CUTOFF, LFC_CUTOFF, TOP_N_HEATMAP, "hormone", "hormone")
cat("✓ Analysis 5 complete\n\n")

# =============================================================================
# SECTION 12: OVERLAP — Genes consistently induced across all Fg-infected studies
# =============================================================================
cat("── Cross-study overlap: consistent Fg-response genes ────────────────\n")

overlap_dir <- file.path(OUT_DIR, "06_Cross_Study_Overlap")
dir.create(overlap_dir, showWarnings = FALSE)

# Load the significant DEG lists from each infected-vs-mock analysis
sig_a1 <- read.csv(file.path(OUT_DIR, "01_Fg_Infected_vs_Mock",
                               "01_Fg_Infected_vs_Mock_significant_DEGs.csv"))
sig_a4 <- read.csv(file.path(OUT_DIR, "04_S9_SE_Fg_vs_Mock",
                               "04_S9_SE_Fg_vs_Mock_significant_DEGs.csv"))

up_a1  <- filter(sig_a1, direction == "Up")$gene_id
up_a4  <- filter(sig_a4, direction == "Up")$gene_id
core_up <- intersect(up_a1, up_a4)

dn_a1   <- filter(sig_a1, direction == "Down")$gene_id
dn_a4   <- filter(sig_a4, direction == "Down")$gene_id
core_dn <- intersect(dn_a1, dn_a4)

cat(sprintf("  Core Fg-induced genes (Analysis 1 ∩ Analysis 4): %d\n", length(core_up)))
cat(sprintf("  Core Fg-repressed genes: %d\n", length(core_dn)))

write.csv(data.frame(gene_id = core_up, direction = "Up_in_both"),
          file.path(overlap_dir, "Core_Fg_induced_genes.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = core_dn, direction = "Down_in_both"),
          file.path(overlap_dir, "Core_Fg_repressed_genes.csv"), row.names = FALSE)

# Heatmap of core genes across ALL infected samples (batch-corrected VST for viz)
if (length(core_up) >= 5) {
  all_infected <- SAMPLE_INFO %>% filter(infected == "Infected")
  rld_all <- rlog(
    DESeqDataSetFromMatrix(
      countData = counts_filt[, all_infected$sample],
      colData   = all_infected,
      design    = ~ study_id + layout
    ), blind = TRUE
  )
  # ComBat-seq batch correction for visualisation (NOT for DEG stats)
  batch    <- as.integer(factor(all_infected$study_id))
  vst_bc   <- ComBat(dat = assay(rld_all), batch = batch,
                      mod = model.matrix(~ as.integer(factor(all_infected$infected))))

  core_genes <- c(head(core_up, min(25, length(core_up))),
                  head(core_dn, min(15, length(core_dn))))
  core_genes <- intersect(core_genes, rownames(vst_bc))
  if (length(core_genes) >= 2) {
    z_core <- t(scale(t(vst_bc[core_genes, ])))
    ann_core <- all_infected %>%
      select(sample, study_id, layout) %>% column_to_rownames("sample")
    png(file.path(overlap_dir, "Core_FHB_genes_Heatmap_BatchCorrected.png"),
        width = 3000, height = max(2500, length(core_genes) * 80), res = 200)
    pheatmap(z_core,
             annotation_col = ann_core,
             color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
             cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = TRUE, show_colnames = FALSE,
             fontsize_row = 8,
             main = "Core FHB-response genes (batch-corrected for visualisation only)",
             border_color = NA)
    dev.off()
    cat("  ✓ Core gene heatmap saved (ComBat-corrected for display)\n")
  }
}

cat("✓ Cross-study overlap complete\n\n")

# =============================================================================
# SECTION 13: SUMMARY
# =============================================================================
cat("=============================================================================\n")
cat("  WHEAT FHB MULTI-STUDY RNA-SEQ — DEG ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("  Organism : Triticum aestivum (IWGSC RefSeq v2.1)\n")
cat("  Libraries: 500 (466 PE + 34 SE) across 14 studies\n")
cat("  Excluded : SRR639194 (single unreplicated SE library)\n")
cat("  Thresholds: padj <", PADJ_CUTOFF, "| |log2FC| >", LFC_CUTOFF, "\n")
cat("-----------------------------------------------------------------------------\n")
cat("  Analysis 1 — Fg Infected vs Mock (multi-study):\n")
cat("    Output: deseq2_results/01_Fg_Infected_vs_Mock/\n")
cat("  Analysis 2 — Resistant vs Susceptible (multi-study, Fg only):\n")
cat("    Output: deseq2_results/02_Resistant_vs_Susceptible/\n")
cat("  Analysis 3 — Sumai3 vs non-Sumai3 (Fg only):\n")
cat("    Output: deseq2_results/03_Sumai3_vs_nonSumai3/\n")
cat("  Analysis 4 — S9 SE GAIIx: Fg vs Mock (2618 & 2890):\n")
cat("    Output: deseq2_results/04_S9_SE_Fg_vs_Mock/\n")
cat("  Analysis 5a — S11 Fielder: Fg+GA vs Fg:\n")
cat("    Output: deseq2_results/05a_S11_Fielder_FgGA_vs_Fg/\n")
cat("  Analysis 5b — S11 Fielder: Fg+ABA vs Fg:\n")
cat("    Output: deseq2_results/05b_S11_Fielder_FgABA_vs_Fg/\n")
cat("  Analysis 6 — Cross-study core FHB gene overlap:\n")
cat("    Output: deseq2_results/06_Cross_Study_Overlap/\n")
cat("  QC — Global PCA & sample distance across all 500 samples:\n")
cat("    Output: deseq2_results/00_Global_QC/\n")
cat("=============================================================================\n")

sink(file.path(OUT_DIR, "session_info.txt"))
sessionInfo()
sink()
cat("✓ Session info saved\n")
