# WADE → Nextflow Conversion Spec Sheet

**Target audience:** an engineer or coding model porting the [`phac-nml/wade`](https://github.com/phac-nml/wade) R package into a Nextflow (DSL2 / nf-core-style) pipeline.

**Source analyzed:** `~/wade` — WGS Analysis and Detection of Molecular Markers (WADE), National Microbiology Laboratory, Public Health Agency of Canada. Apache-2.0.

**Goal:** reproduce WADE's molecular-typing outputs (LabWare upload CSVs + extracted gene FASTAs) as a portable, parallel, reproducible Nextflow pipeline that runs per-sample instead of looping inside R, with no Shiny UI and no hard-coded Windows paths.

---

## Table of contents

1. [What WADE does (domain)](#1-what-wade-does)
2. [Source architecture & control flow](#2-source-architecture)
3. [The organism × analysis matrix](#3-org-analysis-matrix)
4. [Algorithm decomposition (the units that become processes)](#4-algorithm-decomposition)
5. [Inputs](#5-inputs)
6. [Reference data layout (`wade-data`)](#6-reference-data-layout)
7. [Outputs](#7-outputs)
8. [Target Nextflow architecture](#8-target-nextflow-architecture)
9. [The R refactor contract (critical)](#9-the-r-refactor-contract)
10. [Containerization](#10-containerization)
11. [Config, params & schema](#11-config-params-schema)
12. [Behavior-preservation notes & known quirks](#12-behavior-preservation-notes)
13. [Migration plan (phased)](#13-migration-plan)
14. [How to test the result](#14-how-to-test-the-result)
15. [Proposed repo layout](#15-proposed-repo-layout)

---

## 1. What WADE does

WADE takes **pre-assembled bacterial genome contigs** (FASTA) and uses **NCBI BLAST+** to detect and characterize genes for molecular typing of four organisms: **GAS** (Group A *Streptococcus*), **GBS** (Group B *Streptococcus*), **PNEUMO** (*Streptococcus pneumoniae*), **GONO** (*Neisseria gonorrhoeae*).

Per locus it determines: presence/absence, the matching **allele number** (vs a curated multi-FASTA lookup), **amino-acid mutations/motifs** (vs a wild-type reference + a mutations list), whether a gene is **disrupted/pseudogene**, and assembles a per-sample **molecular profile**. Higher-level analyses combine loci into **MLST/NGSTAR/NGMAST sequence types**, **emm types**, **serotypes**, **23S rRNA mutated-allele counts** (from VCF), and finally **predicted MICs / AMR interpretations** formatted for **LabWare** LIMS upload.

The original entry point is an **R Shiny app** (`WADE.R`); a non-interactive batch entry exists in `LinuxWADE.R`. The pipeline logic lives in an installable R package (`R/*.R`, exported via `NAMESPACE`).

---

## 2. Source architecture

```
~/wade
├── WADE.R                # Shiny UI + server: dispatches Org×Test → *_pipeline() calls  (REPLACE)
├── LinuxWADE.R           # Non-interactive batch example (closest to a CLI)             (REFERENCE)
├── DESCRIPTION/NAMESPACE # R package metadata; 48 exported functions
├── R/                    # The actual logic (port this):
│   ├── CommonFunctions.R #   BLAST engine + parsing + LabWare helpers (the core)
│   ├── MASTER.R          #   MASTER_pipeline(): presence/allele/mutation per locus
│   ├── MLST.R            #   MLST_pipeline(): MLST / NGSTAR / NGMAST sequence typing
│   ├── EMM.R / EMM2.R    #   emm typing + M1UK/M1DK variant typing
│   ├── SerotypeR.R       #   SEROTYPE_pipeline(): 4-stage CPS serotyping (PNEUMO/GBS)
│   ├── 16SrRNA.R         #   rRNA16S_pipeline(): species ID from 16S
│   ├── 23SrRNA.R         #   rRNA23S_pipeline(): count mutated 23S alleles FROM VCF (no BLAST)
│   ├── Makeblastdb.R     #   MakeblastDB_pipeline(): index allele lookup FASTAs
│   ├── LabwareUpload_*.R #   GAS/GBS/GONO/PNEUMO AMR/TOXINS/VIRULENCE → LabWare CSV
│   ├── LabwareUpload_Metrics.R # assembly-QC TSV combiner (interactive file picker)
│   └── (Curator utils)   #   UpdateLookups, FindDuplicates, Proximity_Test, MIC_Compare …
├── WADE/                 # Runtime scaffold:
│   ├── DirectoryLocations.csv  # per-organism path config (REPLACE with channels/params)
│   ├── list.csv                # sample list (REPLACE with samplesheet)
│   ├── Example_contigs/        # GAS/GBS/GONO/PNEUMO.fasta + GONO.vcf + PNEUMO.vcf  ← TEST DATA
│   ├── Output/  temp/          # runtime scratch (REPLACE with work dirs / publishDir)
└── wade-data/            # git submodule (phac-nml/wade-data) — reference DBs. NOT cloned here.
```

**Control flow today (per analysis):**
`Shiny/Linux entry` → `*_pipeline(Org, Test, Sample, Locus, curr_work_dir)` → `getdirectory()` reads `DirectoryLocations.csv` → `refdirectory()` builds reference paths → **loop over samples in `list.csv`** → **loop over loci in `loci.csv`** → `system("makeblastdb …")` / `system("blastn …")` → parse text → write `Output/output_profile_<Org>_<Test>.csv` and FASTAs → a `labware_*()` post-processor reads that CSV and writes `LabWareUpload_*.csv`.

**The single most important structural change:** the **sample loop becomes Nextflow's per-sample scatter**, and the **locus loop stays inside the per-sample R call** (or becomes a finer scatter — see §8). Everything else is plumbing.

---

## 3. Org × analysis matrix

Reconstructed from the `switch()` blocks in `WADE.R`. "Routine `ALL`" is what a lab runs per organism; reproduce these as the default per-organism subworkflows.

| Org | Analysis (UI key) | Underlying calls (in order) |
|-----|-------------------|------------------------------|
| **GAS** | AMR | `MASTER(AMR)` → `labware_gas_amr` |
| | EMM | `EMM_pipeline` |
| | M1UK | `EMM_V_pipeline(variant)` |
| | MLST | `MLST(MLST)` |
| | rRNA16S | `MASTER(rRNA16S)` → `rRNA16S_pipeline` |
| | TOXINS | `MASTER(TOXINS)` → `labware_gas_toxins` |
| | **ALL** | AMR+labware, MLST, EMM, rRNA16S+labware, TOXINS+labware |
| **GBS** | AMR | `MASTER(AMR)` → `labware_gbs_amr` |
| | MLST | `MLST(MLST)` |
| | rRNA16S | `MASTER(rRNA16S)` → `rRNA16S_pipeline` |
| | SERO | `SEROTYPE_pipeline` |
| | **ALL** | AMR+labware, MLST, rRNA16S+labware, SERO |
| **PNEUMO** | AMR_ALL | `MASTER(AMR)` → `rRNA23S` → `labware_pneumo_amr` |
| | AMR | `MASTER(AMR)` (alleles only) |
| | rRNA16S | `MASTER(rRNA16S)` → `rRNA16S_pipeline` |
| | MLST | `MLST(MLST)` |
| | rRNA23S | `rRNA23S_pipeline` (VCF) |
| | VIRULENCE | `MASTER(VIRULENCE)` → `labware_pneumo_virulence` |
| | SERO | `SEROTYPE_pipeline` |
| | **ALL** | AMR_ALL+rRNA23S+labware, MLST, VIRULENCE+labware, SERO |
| **GONO** | AMR_ALL | `MASTER(AMR_ALL)` → `MLST(NGSTAR_AMR)` → `rRNA23S` → `labware_gono_amr` |
| | AMR | `MASTER(AMR)` (alleles only) |
| | rRNA23S | `rRNA23S_pipeline` (VCF) |
| | MLST / NGSTAR / NGMAST | `MLST(<scheme>)` |
| | **ALL** | AMR_ALL chain + MLST + NGSTAR + NGMAST |
| **Curator** | UPDATE_LOOKUPS / CONTAMINATION_CHECK / MIC_CHECK / REMOVE_DUPLICATES | maintenance utilities — **out of scope** for v1 pipeline (port later as standalone entry workflows) |
| **all** | LW_METRICS | `metrics()` — assembly-QC combiner; interactive, see §12 |

Note the `Test_id` remaps inside the R (preserve them): `AMR_ALL`→`AMR` (and GONO NGSTAR loci `penA/mtrR/porB/ponA/gyrA/parC/rRNA23S` route to `NGSTAR`); `NGSTAR_AMR`→`NGSTAR` (mutation mode); MakeblastDB maps `AMR_ALL`→`AMR`, `rRNA23S`→`NGSTAR`, and **skips** indexing for `rRNA23S`/`TOXINS`.

---

## 4. Algorithm decomposition

These are the reusable computational units. Each becomes a Nextflow **module** (process). Source: `R/CommonFunctions.R` unless noted.

### 4a. Reference indexing — `quickBlastIndex()` / `MakeblastDB_pipeline()`
For each locus in `loci.csv`, `makeblastdb -in <Lkup_Dir>/<locus>.fasta -dbtype nucl`. One-time per reference set. → process `MAKE_BLAST_DB` (can also be baked into the reference-data container, or run once and shared as a value channel).

### 4b. Per-sample contig indexing — (inline in `MASTER_pipeline`)
`makeblastdb -in <contig> -dbtype nucl`. → fold into each per-sample process, or a dedicated `INDEX_CONTIG` process whose output (the `.n*` files) fans into downstream BLAST processes.

### 4c. Presence/absence + allele + mutation — `MASTER_pipeline()` (the workhorse)
Per sample × per locus:
1. `blasthits()` — `blastn -query <wildgene> -db <contig>` → grep "No hits found" → **POS/NEG**.
2. If POS: `parseblast()` — extract WT vs subject DNA, `Identities`, `%`, coverage.
3. `AAconvert()` — `Biostrings::translate()` WT & query; trim trailing `*`.
4. `pairwiseAlignment(..., substitutionMatrix="BLOSUM50", gapOpening=-2, gapExtension=-8)` → AA alignment → mutation string `S81F`-style.
5. `getmotifs()` — map mutations to named motifs using `reference/loci_mutations.csv` (`Locus_id,Name,Posn_1,Posn_2,WildType`).
6. Disruption call: contains `*` OR coverage < 80 → `"Disrupted"`.
7. Allele lookup: write extracted gene → `blastn` vs `allele_lkup_dna/<locus>.fasta -num_alignments 1` → exact length+100% → allele parts from header `>Locus_AlleleNum_Mutation_Comment` (split on `_`); else `NF`/`???`.
8. Special cases: `porB` coverage <90→`porB1a`, ≥90→`porB1b`; coverage ≤ `wt_id` threshold→force `NEG`.
9. Emit 7 columns/locus: `<locus>_result,_allele,_mutations,_comments,_BlastID,_PctWT,_motifs`, plus a concatenated `SampleProfile` (e.g. `gyrA S81F-parC S79F`, else `Susceptible`).
10. Append extracted sequences to `output_dna.fasta`, `output_dna_notfound.fasta`, `output_aa.fasta`; dedupe not-found → `output_dna_notfound_distinct.fasta`.

→ process `MASTERBLASTR(meta, contig, refdir, test)` emitting one **profile row CSV** + the FASTAs. (Test ∈ AMR/TOXINS/VIRULENCE/rRNA16S/MASTERBLASTR.)

### 4d. Sequence typing — `MLST_pipeline()` + `locusblast()` / `mutationblast()`
Per sample × per locus: `blastquery()` (`blastn … -outfmt 6 -num_alignments 100 -evalue <allele>`) → filter `Ident==100 & Mismatches==0 & Gaps==0` (MLST also `Align > size-10`) → allele number from header. Combine loci → `left_join` against `reference/profiles.txt` to assign **ST**. NGMAST: anti-join `bad_alleles.csv`, prefix `porB-`/`tbpB-`/`ST-`. NGSTAR_AMR: uses `mutationblast()` (mutation column, not ST). → process `MLST(meta, contig, refdir, scheme)`.

### 4e. emm typing — `EMM_pipeline()` / `EMM_V_pipeline()` (M1UK/M1DK SNP variant). → process `EMM`.

### 4f. Serotyping — `SEROTYPE_pipeline()` (PNEUMO/GBS), 4 stages: locus presence → intact/pseudogene → serotype-determining AA substitutions → whole-gene allele match. → process `SEROTYPE`.

### 4g. 23S rRNA from VCF — `rRNA23S_pipeline()` (**no BLAST**)
Reads `<vcf_dir>/<sample>.vcf`, greps organism-specific position lines (GONO `23S4_NCCP11945` 2045/2597; PNEUMO `23S_rRNA_R6_sprr02` 2061/2613), computes alt-allele fraction `AO/DP`, bins into mutated-allele count 0–4 (`<14→0, 15-34→1, 35-64→2, 65-84→3, ≥85→4`). → process `RRNA23S(meta, vcf, org)`. Missing VCF → count `NA` (R) / 0 fallback.

### 4h. 16S species ID — `rRNA16S_pipeline()`: post-processes the `MASTER(rRNA16S)` profile CSV; `<97%`→`NF`. → process `RRNA16S_FORMAT`.

### 4i. LabWare formatters — `labware_gas_amr / gas_toxins / gbs_amr / gono_amr / pneumo_amr / pneumo_virulence`
Consume the aggregated `output_profile_<Org>_<Test>.csv` (+ 23S counts) and emit `LabWareUpload_*.csv`: build molecular profiles via `posneg_gene()` / `allele1SNP()…allele5SNPs()`, compute MICs via `MICcalc()` (`2^calc` rounded; `<=`/`>=`/exact + R/I/DS interpretation). → processes `LABWARE_<ORG>_<TEST>` (run **after** the per-sample collate).

### 4j. Collate — replaces R's `rbind` accumulation: gather per-sample rows into one CSV (`collectFile`/`COLLATE`).

---

## 5. Inputs

| Input | Original source | Nextflow form |
|-------|-----------------|---------------|
| Assembly contigs | `<ContigsDir>/<SampleNo>.fasta` (ext `.fasta`) | samplesheet `assembly` column → `Channel.fromPath` |
| VCFs (GONO/PNEUMO 23S only) | `<VCFDir>/<SampleNo>.vcf` | samplesheet `vcf` column (optional) |
| Sample list | `WADE/list.csv` (`SampleNo,Variable`) | `--input samplesheet.csv` (`sample,assembly,vcf,variable`) |
| Per-org path config | `DirectoryLocations.csv` | **deleted** → `params` + channels |
| Reference DBs | `wade-data/` submodule | `--wade_data <dir>` (staged dir or submodule); see §6 |
| BLAST sensitivity | `reference/blast_evalues.csv` (`contig,allele,wt_id` = `1e-50,1e-98,10`) | keep as a staged ref file or expose as params |
| Which analysis | Shiny radio buttons | `--organism`, `--analysis` (or `all`) |

**Proposed samplesheet** (`assets/samplesheet.csv`):
```csv
sample,assembly,vcf,variable
GAS,/path/Example_contigs/GAS.fasta,,
PNEUMO,/path/Example_contigs/PNEUMO.fasta,/path/Example_contigs/PNEUMO.vcf,
GONO,/path/Example_contigs/GONO.fasta,/path/Example_contigs/GONO.vcf,
```
`variable` is a free-text passthrough (e.g. a phenotypic MIC) carried into outputs via the meta map.

---

## 6. Reference data layout

`wade-data` is a **git submodule** (`https://github.com/phac-nml/wade-data.git`) and is **not cloned** in `~/wade`. The porting model must obtain it (or stage equivalents). Path conventions are derived from `getdirectory()`/`refdirectory()`:

```
<wade_data>/<ORG>/<TEST>/
├── allele_lkup_dna/        # <locus>.fasta  multi-FASTA allele DBs (optional per locus)
│                           #   header: >Locusname_AlleleNum_Mutation_Comment_Variable
│                           #   e.g.    >folA_28_I100L_Sample5_NA
├── wildgenes/              # <locus>.fasta  wild-type reference gene (REQUIRED for MASTER)
└── reference/
    ├── loci.csv            # list of loci (+ size/min cols for MLST/NGMAST)
    ├── profiles.txt        # tab-sep ST→allele profiles (MLST/NGSTAR/NGMAST)
    ├── blast_evalues.csv   # contig,allele,wt_id
    ├── loci_mutations.csv  # optional: Locus_id,Name,Posn_1,Posn_2,WildType
    └── bad_alleles.csv     # NGMAST only
```
`<ORG>` ∈ {GAS, GBS, GONO, PNEUMO}; `<TEST>` ∈ {AMR, TOXINS, VIRULENCE, MLST, NGSTAR, NGMAST, rRNA16S, SERO, …}. Stage this as a single channel (a directory input per `(org,test)`), so each analysis process receives exactly its `reference/` + `wildgenes/` + `allele_lkup_dna/` subtree. **Pin the submodule commit** for reproducibility (the repo auto-updates it daily via `.github/workflows/update-wade-data.yml` — do not track a moving target in releases).

---

## 7. Outputs

Reproduce these (names are LabWare contracts — keep them). Route via `publishDir`/`outdir`.

| File | Producer | Notes |
|------|----------|-------|
| `output_profile_<Org>_<Test>.csv` | MASTER/MLST | per-locus 7-col block + `SampleProfile`. Intermediate. |
| `LabWareUpload_<Org>_<Test>.csv` | labware_* | **primary deliverable** per analysis |
| `LabWareUpload_<Org>_<Test>_good.csv` / `_bad.csv` | MLST/NGMAST, metrics | split by ST-found / QC |
| `output_dna.fasta` | MASTER | all extracted gene DNA |
| `output_aa.fasta` | MASTER | translated proteins |
| `output_dna_notfound[_distinct].fasta` | MASTER | novel alleles (allele=`NF`) |
| `LabWareUpload_METRICS[_good/_bad].csv` | metrics | assembly QC |

Add modern niceties not in the original: a MultiQC-style summary, `pipeline_info/` (trace/timeline/report), and software versions (`versions.yml`).

---

## 8. Target Nextflow architecture

**Parallelization model:** scatter on **sample** (each contig is independent → free parallelism + `-resume`). Within a sample, the **per-locus loop stays inside the R process** for v1 (simplest, preserves R logic). *Optional v2:* scatter per `(sample, locus)` for very large locus sets, then gather — only if profiling shows it's worth the extra join complexity. Recommend **v1 first**.

**Dataflow (per organism `ALL` routine), e.g. GAS:**
```
samplesheet ─► parse ─► ch_samples = [ meta(id,org,variable), assembly, vcf? ]

ch_samples ┬─► MASTERBLASTR(AMR)    ─► per-sample row ─┐
           ├─► MLST(MLST)           ─► per-sample row ─┤
           ├─► EMM                  ─► per-sample row ─┤
           ├─► MASTERBLASTR(rRNA16S)─► RRNA16S_FORMAT ─┤
           └─► MASTERBLASTR(TOXINS) ─► per-sample row ─┤
                                                       ▼
                            COLLATE (collectFile by org,test)
                                                       ▼
                  LABWARE_GAS_AMR / _TOXINS / RRNA16S formatting
                                                       ▼
                                   publishDir(outdir)
```

**Module inventory (`modules/local/`):**
- `makeblastdb.nf` — index reference allele DBs (or contig)
- `masterblastr.nf` — §4c; param `test`
- `mlst.nf` — §4d; param `scheme` (MLST/NGSTAR/NGMAST/NGSTAR_AMR)
- `emm.nf`, `emm_variant.nf` — §4e
- `serotype.nf` — §4f
- `rrna23s.nf` — §4g (VCF)
- `rrna16s_format.nf` — §4h
- `labware_gas_amr.nf`, `labware_gas_toxins.nf`, `labware_gbs_amr.nf`, `labware_gono_amr.nf`, `labware_pneumo_amr.nf`, `labware_pneumo_virulence.nf`
- `collate.nf` — gather per-sample rows
- `custom/dumpsoftwareversions.nf`

**Subworkflows (`subworkflows/local/`):** `gas.nf`, `gbs.nf`, `gono.nf`, `pneumo.nf` — each takes `ch_samples` + `ch_refdata`, wires its analyses (mirror §3), `emit:` LabWare channels. Top `workflows/wade.nf` branches by `meta.org` (or `--organism`) using `.branch{}`.

**meta map convention:** `[ id:'GAS-001', org:'GAS', variable:'8ug/ml' ]` carried in every tuple alongside files (nf-core style). VCF optional → use `[meta, assembly, vcf]` with `vcf` possibly `[]`.

---

## 9. The R refactor contract

This is where most of the work is. **Do not rewrite the science — wrap it.** Refactor the package so each analysis runs **non-interactively on one sample** with **explicit paths**, then call it from `bin/`.

**Required changes to the R:**
1. **Remove the UI/host coupling:** drop `library(shiny/shinyWidgets/DT/beepr)`; delete `WADE.R` Shiny app, `shell.exec`, `beep()`, and `choose.files()` (in `metrics()` — replace with CLI args).
2. **Replace `getdirectory()`:** stop reading `DirectoryLocations.csv`. Accept `--contig`, `--vcf`, `--refdir`, `--outdir`, `--tempdir` (default tempdir to the Nextflow task dir / `tempdir()`).
3. **Remove the sample loop:** `MASTER_pipeline`/`MLST_pipeline`/etc. currently `for(m in 1:NumSamples)` over `list.csv`. Make the core operate on a **single** `SampleNo`/contig; Nextflow supplies one sample per task. Keep the **locus loop** inside.
4. **Single-row output:** emit a one-row profile CSV per sample (so `COLLATE` can concatenate). Keep exact column names/order.
5. **POSIX paths:** drop `C:\WADE\` assumptions and `\\` separators; use `file.path()`.
6. **`system("makeblastdb"/"blastn")` stays** — BLAST+ comes from the container; make sure the contig DB is built in the task workdir.
7. **Provide a CLI dispatcher** `bin/wade.R` (use `optparse` or `argparser`), e.g.:
   ```
   wade.R masterblastr --org GAS --test AMR --contig in.fasta --refdir REF/ --outdir .
   wade.R mlst         --org GONO --scheme NGSTAR --contig in.fasta --refdir REF/ --outdir .
   wade.R rrna23s      --org PNEUMO --vcf in.vcf --outdir .
   wade.R labware      --org GAS --test AMR --profile output_profile_GAS_AMR.csv --outdir .
   ```
   Each subcommand sources the (refactored) package functions. Mark `bin/*.R` executable with a `#!/usr/bin/env Rscript` shebang so Nextflow finds them on `PATH`.
8. Keep `CommonFunctions.R` helpers (`parseblast`, `AAconvert`, `getmotifs`, `locusblast`, `mutationblast`, `posneg_gene`, `allele{1..5}SNP`, `MICcalc`) **byte-for-byte where possible** — they encode validated lab logic.

Keep the R package installable (it already passes `R CMD check` in CI) so functions stay unit-testable in R independently of Nextflow.

---

## 10. Containerization

One image with **R + Bioconductor + tidyverse + BLAST+**. Biostrings/tidysq are the only non-trivial deps.

- **Base:** `bioconductor/bioconductor_docker:RELEASE_3_19` (or `rocker/r-ver` + `BiocManager`).
- **R pkgs:** `plyr, tidyverse, tidyselect, tidysq, data.table, optparse` (+ `Biostrings` from Bioconductor). Drop `shiny/shinyWidgets/DT/beepr/readxl` unless `metrics()` still needs `readxl`.
- **System:** `ncbi-blast+` (`apt-get install -y ncbi-blast+`, same as `.github/workflows/ci.yml`).
- Provide **both** a `Dockerfile` and a Conda `environment.yml` (`bioconda::blast`, `bioconda::bioconductor-biostrings`, `conda-forge::r-tidyverse`, …) so `-profile docker|singularity|conda` all work. Consider **Wave** to build from the conda env automatically.
- Pin versions (R, Bioconductor release, BLAST+) and emit them to `versions.yml`.
- Set `ENV BLASTDB_LMDB_MAP_SIZE=1000000` (the documented Windows workaround; harmless on Linux).

---

## 11. Config, params & schema

- `nextflow.config`: `params { input; outdir; organism; analysis='all'; wade_data; blast_evalue_contig='1e-50'; blast_evalue_allele='1e-98'; wt_id_threshold=10 }`; `manifest`; `profiles { test; docker; singularity; conda }`.
- `conf/base.config`: resource labels (`process_single/low/medium`) + `errorStrategy 'retry'` with `task.attempt` scaling. BLAST per locus is light; right-size small.
- `conf/test.config`: points `input` at a samplesheet built from `WADE/Example_contigs/` and a pinned `wade_data` test subset.
- `conf/modules.config`: `publishDir` + `ext.args` for BLAST flags (don't hard-code evalues in the process script — pass via `ext.args`).
- `nextflow_schema.json`: build with `nf-core pipelines schema build`; validate with `nf-validation`/`nf-schema`.

---

## 12. Behavior-preservation notes

Port these **exactly**; they look like bugs or oddities but define current (validated) output. Add a regression test for each.

1. **`locusblast()` MLST size filter fall-through** (`CommonFunctions.R` ~L256-268): the `if (MLST) {…} else if (NGMAST) {…}` is immediately followed by a bare `{ df.blastout100 <- filter(Ident==100 & Mismatches==0 & Gaps==0) }` block that **always re-runs**, overriding the `Align > size-10` filter. Preserve observed behavior; do **not** "fix" unless the lab confirms intent.
2. **Allele header parsing** depends on `_`-delimited fields (`AlleleParts[2]`/`[3]`/`[4]`; NGMAST uses `[3]`). Reference FASTA headers must keep the `Locus_Allele_Mutation_Comment_Variable` convention.
3. **Disruption rule:** stop codon `*` in AA alignment **or** `IDcoverage < 80`.
4. **`wt_id` threshold** (default 10) forces `NEG` when `IDcoverage <= threshold`.
5. **porB special-casing** (`porB1a` <90%, `porB1b` ≥90%).
6. **AA translate** trims a single trailing `*`; alignment uses **BLOSUM50, gapOpen −2, gapExtend −8** — keep exact.
7. **23S VCF bins** (0–4) and the **exact grep strings** per organism — keep literally; they match specific VCF reference headers.
8. **`metrics()` is interactive** (`choose.files`) and combines two QC TSVs (`Assembly_stats.tsv`, `fastqc.tsv`) — refactor to take two `--file` args; thresholds: N50<30k→Warning, Coverage<10 | N50<10k | MeanContig<1000→Fail.
9. **CSV quirks:** outputs are written `quote = FALSE` — match quoting/whitespace or your golden-file diff will flag spurious mismatches (normalize before comparing, see §14).
10. **Curator utilities** (UpdateLookups/Proximity/MIC_Compare/FindDuplicates) are admin tools, not per-sample analyses — defer to a later phase as separate entry workflows.

---

## 13. Migration plan (phased)

1. **Phase 0 — scaffold:** `nf-core pipelines create`; add container; obtain & pin `wade-data`; build samplesheet from `WADE/Example_contigs/`.
2. **Phase 1 — R refactor:** make `MASTERBLASTR` single-sample + CLI (`bin/wade.R`). Get one process running on one contig; compare to R output. *This de-risks everything.*
3. **Phase 2 — MASTER path end-to-end:** `INDEX → MASTERBLASTR(AMR) → COLLATE → LABWARE_GAS_AMR`. Validate vs golden file.
4. **Phase 3 — remaining analyses:** MLST/NGSTAR/NGMAST, EMM, SERO, rRNA16S, rRNA23S(VCF) + their LabWare formatters.
5. **Phase 4 — subworkflows:** wire per-org `ALL` routines (§3); branch by org.
6. **Phase 5 — hardening:** nf-test suite, `-profile test`, CI, lint, schema, docs, `-resume`/stub correctness, MultiQC/versions.
7. **Phase 6 (optional):** Curator utilities; per-locus scatter (v2).

---

## 14. How to test the result

A typing pipeline is only useful if outputs are **scientifically identical** to the validated R tool. Use five layers, cheap → strict:

### Layer 1 — Plumbing / stub (`-stub-run`)
Give every process a `stub:` block that `touch`es expected outputs. Validates channel wiring/joins with no BLAST/R:
```bash
nextflow run . -profile test --stub-run --outdir results_stub
```

### Layer 2 — Environment smoke test (bundled `test` profile)
Uses the in-repo example contigs (`WADE/Example_contigs/{GAS,GBS,GONO,PNEUMO}.fasta` + `GONO.vcf`, `PNEUMO.vcf`). Proves containers + BLAST + R work end-to-end on tiny data:
```bash
nextflow run . -profile test,docker   --outdir results       # or singularity / conda
```

### Layer 3 — `nf-test` unit/module/pipeline tests (the regression net)
- **Per module:** feed one example contig, snapshot the output CSV/FASTA (`assert snapshot(...).match()`). Tag by analysis.
- **Subworkflow:** GAS/GBS/GONO/PNEUMO `ALL` produce expected file sets.
- **Pipeline:** full `-profile test` run; assert all `LabWareUpload_*.csv` exist and match snapshots/md5.
```bash
nf-test test --profile docker            # all
nf-test test tests/modules/masterblastr.nf.test
```
Use `nf-core pipelines lint` + `prettier`/`editorconfig` to keep nf-core-compliant.

### Layer 4 — Golden-file parity vs original R (the scientific gate)
The decisive test. Run the **original** WADE on the example contigs, capture its outputs as fixtures, and assert the Nextflow pipeline reproduces them.

1. Generate reference outputs once (Linux path, BLAST+ installed):
   ```bash
   # in ~/wade after staging wade-data + DirectoryLocations.csv for Example_contigs
   Rscript LinuxWADE.R          # writes WADE/Output/output_profile_*.csv, LabWareUpload_*.csv
   ```
   Store under `tests/fixtures/expected/`.
2. Run the Nextflow pipeline on the same inputs → `results/`.
3. **Normalize then diff** (outputs are unquoted, column order/whitespace/row order can differ). Compare as data, not bytes:
   ```bash
   # sort rows + columns, trim whitespace, then diff — e.g. a tiny csvkit/pandas/R harness
   python tests/bin/compare_csv.py results/LabWareUpload_GAS_AMR.csv \
                                    tests/fixtures/expected/LabWareUpload_GAS_AMR.csv
   ```
   Provide `tests/bin/compare_csv.py` (or an R `waldo::compare`) that: reads both, sorts by `SampleNo` + sorts columns, coerces types, and asserts equality with a clear per-cell diff. Treat **any** mismatch as a failure until a lab scientist signs off (it usually means a §12 quirk wasn't preserved).
4. Wire this as an nf-test that runs after the pipeline, or a CI step.

### Layer 5 — Determinism & cache
- **Reproducibility:** run twice, assert identical outputs (BLAST ties / unsorted reads can reorder rows — that's why Layer 4 sorts).
- **`-resume`:** re-run and confirm tasks are `cached`, outputs unchanged.
- **Negative/edge inputs:** a deliberately bad/empty contig → `Sample_Err` propagates (don't crash the run; use `errorStrategy`); a sample with a missing VCF → 23S count `NA`/0 path.

### CI (GitHub Actions)
Mirror `~/wade/.github/workflows/ci.yml` but for Nextflow:
```yaml
# .github/workflows/ci.yml
jobs:
  test:
    steps:
      - uses: actions/checkout@v4
        with: { submodules: true }          # wade-data
      - uses: nf-core/setup-nextflow@v2
      - run: nextflow run . -profile test,docker --outdir results
      - uses: nf-core/setup-nf-test@v1       # or install nf-test
      - run: nf-test test --profile docker
      - run: nf-core pipelines lint
```
Add the Layer-4 parity job (build the R container once, generate fixtures or restore cached fixtures, diff).

**Minimum bar to call the port "done":** Layers 1–3 green in CI **and** Layer 4 parity passes for every organism's `ALL` routine on the bundled example contigs.

---

## 15. Proposed repo layout

```
wade-nf/
├── main.nf
├── nextflow.config
├── nextflow_schema.json
├── conf/{base,test,modules}.config
├── assets/{samplesheet.csv, schema_input.json}
├── bin/                      # refactored single-sample R CLI (chmod +x, Rscript shebang)
│   ├── wade.R                #   dispatcher
│   └── (analysis scripts or one dispatcher sourcing the WADE pkg)
├── modules/local/*.nf        # §8 inventory
├── subworkflows/local/{gas,gbs,gono,pneumo}.nf
├── workflows/wade.nf
├── containers/{Dockerfile, environment.yml}
├── tests/
│   ├── modules/*.nf.test
│   ├── pipeline.nf.test
│   ├── fixtures/expected/    # golden LabWare CSVs from original R
│   └── bin/compare_csv.py
├── wade-data/                # submodule (pinned) OR --wade_data param
└── .github/workflows/ci.yml
```

---

### Quick-reference: source files the porting model must read first
1. `R/CommonFunctions.R` — the BLAST engine + all helpers (port verbatim where possible).
2. `R/MASTER.R` — the per-locus presence/allele/mutation loop.
3. `R/MLST.R` — sequence typing + ST lookup.
4. `WADE.R` (server `switch` blocks) — the authoritative Org×Test routing (§3).
5. `LinuxWADE.R` — the non-interactive call pattern to mimic in `bin/`.
6. `R/23SrRNA.R` — the only VCF-based, non-BLAST analysis.
7. `standalones/*/*_Readme.md` — plain-English algorithm descriptions of MasterBlastR / SerotypeR / WamR-Pneumo.
