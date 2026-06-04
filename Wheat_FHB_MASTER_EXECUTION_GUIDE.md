# WHEAT FHB MULTI-STUDY RNA-seq — MASTER EXECUTION GUIDE
## Organism: Triticum aestivum | Genome: IWGSC RefSeq v2.1
## Libraries: 501 total (466 PAIRED-end + 35 SINGLE-end) | 15 studies
## Theme: Fusarium Head Blight (FHB) resistance
## Pipeline: nf-core/fetchngs → nf-core/rnaseq v3.14.0 → DESeq2

---

## STUDY INVENTORY (all 501 libraries)

| ID | n | Layout | Instrument | Tissue | Cultivar(s) | Treatment | Key Variable |
|----|---|--------|------------|--------|-------------|-----------|--------------|
| S1 | 12 | PE | HiSeq 2500 | Spikelets | norm | Fg 0/24/48/72h | Time-course |
| S2 | 22 | PE | HiSeq 2000 | Spikelet + Rachis | 260-1-1-2 (Fhb1+) / 260-1-1-4 (Fhb1-) | Fg inoculated | Fhb1 locus effect |
| S3 | 2  | PE | HiSeq 4000 | Rachis | Sumai3 / Pasteur | Fg inoculated | Sumai3 vs susceptible |
| S4 | 12 | PE | HiSeq 2000 | Leaf | Sumai3 | Fg isolates R vs S | Pathogen virulence |
| S5 | 100 | PE | HiSeq 2500 | Whole head | 100 doubled haploid lines | Fg inoculated | QTL discovery |
| S6 | 12 | PE | HiSeq 2500 | Spike | Chinese Spring / CS+7EL | Fg vs mock | Chromosome 7EL |
| **S7** | **1** | **SE** | HiSeq 2000 | — | — | — | **⚠ EXCLUDED — no replicates** |
| S8 | 12 | PE | HiSeq 2500 | Rachis | Chinese Spring | unknown | Reference rachis |
| S9 | 22 | **SE** | GA IIx | Spikelet + Rachis | 2618 / 2890 | Fg vs mock | Cultivar × tissue |
| S10 | 2 | PE | HiSeq 4000 | Spike | Thinopyrum substitution lines | none | Alien chromosome |
| S11 | 15 | PE | HiSeq 2500 | Spikelet | Fielder | Fg / Fg+GA / Fg+ABA | Hormone interaction |
| S11b | 6 | PE | HiSeq 2000 | Spikelet | Xiaoyan 22 | Fg osp24 mutant vs PH1 wt | Fg mutant virulence |
| S12 | 192 | PE | HiSeq 2500 | Spike + Rachis | 96 European cultivars × 2 reps | Fg inoculated | FHB resistance panel |
| S13 | 12 | **SE** | BGISEQ-500 | Spikelet + Rachis | unknown | none | Flowering-stage atlas |
| S14 | 7 | PE | HiSeq X Ten | Wheat grains | Sumai3 / Annong0711 / Zhongmai66B | Fg inoculated | Grain FHB |
| S15 | 72 | PE | HiSeq 2000 | Spikelet | NIL38 / NIL51 | Fg vs mock | NIL FHB QTL |

**S12 FHB resistance group breakdown (192 cultivars, all Fg-inoculated, anthesis):**

| Group | n | Sumai3 lineage |
|-------|---|----------------|
| highly_resistant_Sumai3 | 18 | Sumai3 |
| Resistant | 36 | non-Sumai3 |
| moderate_resistance | 90 | non-Sumai3 |
| susceptible | 36 | non-Sumai3 |
| excluded | 12 | Sumai3 (introgressed lines, excluded from resistance analysis) |

**SE samples (35 total):**
- S9: 22 × GAIIx SE (ERR2275046–ERR2275067) — cultivars 2618 & 2890
- S13: 12 × BGISEQ-500 SE (SRR17023792–SRR17023803) — spikelet & rachis
- S7: 1 × HiSeq 2000 SE (SRR639194) — **EXCLUDED**

**Batch correction approach:**
- All multi-study DESeq2 models include `study_id` as a covariate
- `library_layout` (PE/SE) included as covariate in cross-study analyses
- SE and PE are balanced across conditions (not confounded) — safe to combine
- ComBat-seq applied to VST matrix **for visualisation only** (heatmaps, PCA)
- Raw counts fed to DESeq2 without prior correction

---

## FILES PROVIDED

| File | Purpose |
|------|---------|
| `Wheat_FHB_ids.csv` | 501 SRR/ERR accessions for nf-core/fetchngs |
| `Wheat_FHB_nextflow.config` | nf-core/rnaseq v3.14.0 pipeline configuration |
| `Wheat_FHB_lifecycle_policy.json` | S3 Glacier archival policy |
| `Wheat_FHB_DESeq2_analysis.R` | Complete DESeq2 analysis (6 analyses + global QC) |
| `Wheat_FHB_MASTER_EXECUTION_GUIDE.md` | This file |
| `Wheat_FHB_QUICK_REFERENCE.md` | All commands + S12 cultivar table |

---

## PHASE 1: EC2 INSTANCE SELECTION

Hexaploid wheat genome (IWGSC RefSeq v2.1) is ~15 GB uncompressed.
STAR genome index requires ~240 GB RAM to build.

**Recommended instances:**

| Phase | Instance | vCPU | RAM | Cost/hr |
|-------|----------|------|-----|---------|
| STAR genome build (one-time) | r5.24xlarge | 96 | 768 GB | ~$6 |
| Alignment (501 samples) | r5.8xlarge | 32 | 256 GB | ~$2 |
| fetchngs download | m5.4xlarge | 16 | 64 GB | ~$0.77 |
| DESeq2 (local RStudio/R) | Any | — | ≥32 GB | — |

**Tip:** Build the STAR index on r5.24xlarge, save it to S3, then use r5.8xlarge for alignment with `-resume`.

---

## PHASE 2: AWS SETUP

```bash
# Apply S3 lifecycle policy
aws s3api put-bucket-lifecycle-configuration \
    --bucket my-bingya-bucket \
    --lifecycle-configuration file://Wheat_FHB_lifecycle_policy.json

# Verify
aws s3api get-bucket-lifecycle-configuration --bucket my-bingya-bucket \
    | grep -A3 "Wheat_FHB"

# Create project directory on EC2
mkdir -p ~/bingya/Wheat_FHB && cd ~/bingya/Wheat_FHB
cp /path/to/Wheat_FHB_ids.csv .
cp /path/to/Wheat_FHB_nextflow.config .
```

---

## PHASE 3: DOWNLOAD ALL 501 LIBRARIES

```bash
cd ~/bingya/Wheat_FHB

nextflow run nf-core/fetchngs -r 1.12.0 \
    --input Wheat_FHB_ids.csv \
    --outdir ./raw_data \
    --nf_core_pipeline rnaseq \
    -profile docker

# Verify: should show 501 rows (1 header + 500 data + SRR639194 which we exclude)
wc -l ./raw_data/samplesheet/samplesheet.csv

# Check SE vs PE auto-detection
awk -F',' 'NR>1{print $4}' ./raw_data/samplesheet/samplesheet.csv | sort | uniq -c
# Expected: ~466 PE (fastq_2 filled) + 35 SE (fastq_2 empty)
```

**Expected download size:** ~1.5–2 TB total (501 libraries × ~3 GB average)
**Estimated time:** 4–8 hours on m5.4xlarge with 10 Gbps network

---

## PHASE 4: DOWNLOAD IWGSC REFERENCE GENOME

```bash
cd ~/bingya/Wheat_FHB

# IWGSC RefSeq v2.1 from Ensembl Plants release 57
# ⚠ FASTA is ~14 GB compressed — ensure 100 GB free disk space
wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz

wget -c https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gtf/triticum_aestivum/Triticum_aestivum.IWGSC.57.gtf.gz

# Verify checksums
md5sum Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
md5sum Triticum_aestivum.IWGSC.57.gtf.gz

# Upload to S3
aws s3 cp Triticum_aestivum.IWGSC.dna.toplevel.fa.gz \
    s3://my-bingya-bucket/Wheat_FHB/genome/ --no-progress
aws s3 cp Triticum_aestivum.IWGSC.57.gtf.gz \
    s3://my-bingya-bucket/Wheat_FHB/genome/ --no-progress

# Verify
aws s3 ls s3://my-bingya-bucket/Wheat_FHB/genome/
```

---

## PHASE 5: RUN nf-core/rnaseq PIPELINE

```bash
cd ~/bingya/Wheat_FHB

# Recommended: run inside a screen session so it survives disconnection
screen -S wheat_rnaseq

nextflow run nf-core/rnaseq -r 3.14.0 \
    -profile docker \
    -c Wheat_FHB_nextflow.config \
    -resume

# Detach: Ctrl+A then D
# Reattach: screen -r wheat_rnaseq
```

**Monitor progress (second terminal):**
```bash
# Live process status
watch -n 60 'nextflow log -f name,status,duration,exit | tail -30'

# EC2 resource usage
watch -n 10 'free -h && df -h /home && echo "---" && top -bn1 | head -15'

# S3 results accumulating
aws s3 ls s3://my-bingya-bucket/Wheat_FHB/results/ --recursive | wc -l
```

**Estimated runtime:** 24–48 hours for 501 samples on r5.8xlarge

**Key output files after completion:**
```
s3://my-bingya-bucket/Wheat_FHB/results/
├── star_salmon/
│   ├── salmon.merged.gene_counts.tsv      ← DESeq2 input
│   ├── salmon.merged.gene_tpm.tsv         ← TPM values
│   └── *.markdup.sorted.bam               ← alignment BAMs
├── multiqc/
│   └── multiqc_report.html                ← QC report (review FIRST)
├── fastqc/                                ← per-sample FastQC
├── trim_galore/                           ← trimming reports
└── pipeline_info/
    └── execution_report.html
```

---

## PHASE 6: REVIEW MultiQC REPORT

**Before running DESeq2**, download and review the QC report:

```bash
aws s3 cp s3://my-bingya-bucket/Wheat_FHB/results/multiqc/multiqc_report.html ./

# Open in browser (from WSL2):
explorer.exe multiqc_report.html
```

**Key metrics to check for 501 mixed-library samples:**

| Metric | Expected range | Action if outside range |
|--------|---------------|------------------------|
| % Uniquely Mapped | ≥ 70% | Flag and consider excluding |
| SE samples (S9/S13) mapping | May be 5–10% lower than PE | Normal — include `layout` covariate |
| Duplication rate | < 60% | >80% may indicate low-complexity library |
| % rRNA reads | < 10% | High rRNA = poor depletion — flag sample |
| GC content | 44–52% for wheat | Outliers may have contamination |
| M reads (depth) | ≥ 10M | Flag samples with < 5M reads |

**Samples to consider flagging/excluding:**
- SRR639194 — already excluded (single SE library, no metadata)
- Any sample with < 60% unique mapping rate
- Any sample with < 5M mapped reads

---

## PHASE 7: DOWNLOAD COUNTS AND PREPARE FOR DESeq2

```bash
mkdir -p ~/bingya/Wheat_FHB/deseq2_input
cd ~/bingya/Wheat_FHB/deseq2_input

# Download count matrix and QC report
aws s3 cp s3://my-bingya-bucket/Wheat_FHB/results/star_salmon/salmon.merged.gene_counts.tsv ./
aws s3 cp s3://my-bingya-bucket/Wheat_FHB/results/multiqc/multiqc_report.html ./

# Verify count matrix dimensions
head -1 salmon.merged.gene_counts.tsv | tr '\t' '\n' | wc -l   # should be ~502 (1 gene col + 501 samples)
wc -l salmon.merged.gene_counts.tsv                             # should be ~100,000+ gene rows (IWGSC v2.1)

# Copy DESeq2 script
cp /path/to/Wheat_FHB_DESeq2_analysis.R .
```

---

## PHASE 8: INSTALL R PACKAGES (first time only)

```bash
R --no-save << 'EOF'
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cran.r-project.org")
BiocManager::install(c("DESeq2", "apeglm", "sva"))
install.packages(c("tidyverse", "pheatmap", "RColorBrewer",
                   "ggrepel", "patchwork", "scales"),
                 repos = "https://cran.r-project.org")
EOF
```

---

## PHASE 9: RUN DESeq2 ANALYSIS

```bash
cd ~/bingya/Wheat_FHB/deseq2_input

# Run analysis (recommended: screen session, ~30–60 min for 501 samples)
screen -S deseq2
Rscript Wheat_FHB_DESeq2_analysis.R 2>&1 | tee deseq2_run.log

# Monitor progress
tail -f deseq2_run.log
```

**Six analyses run automatically:**

| Analysis | Comparison | Samples used | Design formula | Notes |
|----------|-----------|--------------|----------------|-------|
| 0 | Global QC PCA | All 500 | `~ 1` | Diagnose batch effects |
| 1 | Fg Infected vs Mock | ~150 (S2,S6,S9,S11,S15) | `~ study_id + layout + infected` | Core FHB response |
| 2 | Resistant vs Susceptible | S2 + S12 + S15 | `~ study_id + layout + resistance` | FHB resistance genes |
| 3 | Sumai3 vs non-Sumai3 | S2,S3,S14,S15 | `~ study_id + layout + sumai3` | Sumai3 signature |
| 4 | SE-only: Fg vs Mock | 22 (S9) | `~ cultivar + tissue + infected` | SE-only clean analysis |
| 5a/b | Fg+GA / Fg+ABA vs Fg | 15 (S11) | `~ hormone` | Hormone × FHB |
| 6 | Cross-study core genes | Intersection of A1 ∩ A4 | — | Robust FHB genes |

---

## PHASE 10: UPLOAD RESULTS TO S3

```bash
# Upload all DESeq2 outputs
aws s3 sync deseq2_results/ \
    s3://my-bingya-bucket/Wheat_FHB/deseq2_results/ \
    --exclude "*.RData"

# Verify
aws s3 ls s3://my-bingya-bucket/Wheat_FHB/deseq2_results/ --recursive | wc -l
```

---

## TROUBLESHOOTING

| Problem | Cause | Fix |
|---------|-------|-----|
| STAR genome build OOM | Needs 240 GB RAM | Use r5.24xlarge for genome build only |
| Nextflow runs out of disk | work/ dir grows to 500 GB+ | Add `cleanup = true` (already in config) |
| SE samples have lower alignment | Normal for GAIIx reads | Include `layout` covariate in DESeq2 |
| `lfcShrink` coef not found | DESeq2 coefficient name mismatch | Run `resultsNames(dds)` and update coef string |
| ComBat-seq fails | Batch and mod matrix collinear | Reduce batch levels or skip visualisation correction |
| Very slow fetchngs | 501 accessions, large files | Split ids.csv into batches of 100; run sequentially with `-resume` |
| ERR accessions failing | ENA download throttle | Add `--force_sratools_download` flag to fetchngs |
| salmon.merged.gene_counts.tsv columns are not SRR IDs | fetchngs renames columns to sample names | Update `SAMPLE_INFO$sample` in R script to match |

---

## COST ESTIMATE

| Item | Cost |
|------|------|
| r5.24xlarge (STAR genome build, ~2 hr) | ~$12 |
| r5.8xlarge (alignment, 501 samples, ~36 hr) | ~$72 |
| S3 raw_data (~2 TB, Standard → Glacier Day 30) | $46/mo → $7/mo |
| S3 genome (~15 GB, Standard → Glacier Day 60) | $0.35/mo → $0.05/mo |
| S3 results (~50 GB, Standard) | $1.15/mo |
| S3 DESeq2 results (~1 GB) | $0.02/mo |
| **One-time run cost** | **~$84** |
| **Ongoing monthly storage** | **~$8/month** |
