# wade-nf

A Nextflow (DSL2 / nf-core-style) port of [`phac-nml/wade`](https://github.com/phac-nml/wade) —
**WGS Analysis and Detection of molecular markers** (National Microbiology Laboratory, PHAC).

wade-nf reproduces WADE's molecular-typing outputs (LabWare upload CSVs + extracted gene
FASTAs) as a portable, parallel, reproducible pipeline that runs **per-sample** instead of
looping inside R — no Shiny UI, no hard-coded Windows paths. See [`SPEC.md`](SPEC.md) for the
full conversion specification.

## Status

| Organism | Routine | State |
|----------|---------|-------|
| **GAS** | AMR, TOXINS, MLST, EMM, rRNA16S (`ALL`) | ✅ implemented & **parity-validated** vs original R |
| GBS | AMR, MLST, rRNA16S, SERO | 🚧 modules + R helpers in place; subworkflow + labware pending |
| PNEUMO | AMR_ALL+23S, MLST, VIRULENCE, SERO | 🚧 23S/MLST/MASTERBLASTR ported; pneumo labware pending |
| GONO | AMR_ALL+NGSTAR+23S, MLST/NGMAST | 🚧 MLST/NGSTAR/NGMAST/23S ported; gono labware pending |

The GAS `ALL` routine passes the spec's **minimum bar**: Layers 1–3 green and Layer-4 parity
(byte-identical *as data*) against the original WADE on the bundled example contig.

## Quick start

```bash
# 1. Get reference data (pinned submodule) and the dev/runtime environment
git submodule update --init --recursive          # fetches wade-data
conda env create -f environment.yml              # dev superset (nextflow, nf-core, nf-test, R, BLAST+)
conda activate wade-nf

# 2. Smoke test on the bundled GAS example contig
nextflow run . -profile test,local_env --outdir results          # R + BLAST from the active conda env
#   or, containerised:
nextflow run . -profile test,docker  --outdir results
nextflow run . -profile test,conda   --outdir results            # uses containers/environment.yml

# 3. Real run
nextflow run . \
    --input assets/samplesheet.csv \
    --wade_data wade-data \
    --outdir results \
    -profile docker -resume
```

### Profiles

- `test` — bundled example contig + pinned `wade-data` (tiny data; proves the stack works).
- `docker` / `singularity` / `conda` — container engines (mutually exclusive).
- `local_env` — no containerisation; uses R/BLAST already on `PATH` (e.g. the `wade-nf` conda env). Handy for dev/CI.

## Inputs

Samplesheet (`--input`), comma-separated:

```csv
sample,organism,assembly,vcf,variable
GAS,GAS,GAS.fasta,,
PNEUMO,PNEUMO,PNEUMO.fasta,PNEUMO.vcf,
```

- `organism` ∈ {GAS, GBS, GONO, PNEUMO}; if blank, `--organism` is used. *(This column is a
  practical addition to the spec's 4-column sheet so a single run can mix organisms; the pipeline
  branches on it.)*
- `assembly` — assembled contigs FASTA. Relative paths resolve against the samplesheet's directory.
- `vcf` — optional, only used for 23S rRNA (GONO/PNEUMO).
- `variable` — free-text passthrough carried into outputs via the meta map.

Reference DBs come from `--wade_data` (the `phac-nml/wade-data` layout: `<ORG>/<TEST>/{wildgenes,allele_lkup_dna,reference}`).

## Outputs (`--outdir`)

```
results/
├── GAS/
│   ├── LabWareUpload_GAS_AMR.csv         # primary deliverables (LabWare contracts)
│   ├── LabWareUpload_GAS_TOXINS.csv
│   ├── LabWareUpload_GAS_emm.csv
│   ├── LabWareUpload_GAS_rRNA16S.csv
│   ├── LabWareUpload_GAS_MLST{,_good,_bad}.csv
│   └── per_sample/output_{dna,aa}.fasta  # extracted gene DNA / protein
└── pipeline_info/                        # trace / timeline / report / DAG / software_versions
```

## Architecture

```
samplesheet ─► parse ─► [meta(id,org,variable), assembly, vcf?]
                          │ branch by org
                          ▼
   per-sample scatter ──► MASTERBLASTR(AMR/TOXINS/rRNA16S), MLST, EMM ─► per-sample row CSV
                          │
                          ▼  collectFile (replaces R rbind)
            collated output_profile_<Org>_<Test>.csv
                          │
                          ▼  post-collate formatters
       LABWARE_GAS_AMR / _TOXINS / RRNA16S_FORMAT / MLST_SPLIT ─► publishDir
```

- **`bin/wade.R`** — single-sample, path-explicit CLI dispatcher (`optparse`). Replaces the Shiny
  app / `DirectoryLocations.csv` / `list.csv`. Subcommands: `masterblastr`, `mlst`, `emm`,
  `rrna23s`, `rrna16s`, `labware`, `mlst_split`.
- **`bin/lib/`** — the refactored R. `common.R` keeps the validated helpers (`parseblast`,
  `AAconvert`, `getmotifs`, `locusblast`, `mutationblast`, `posneg_gene`, `allele{1..5}SNP(s)`,
  `MICcalc`) **byte-for-byte**; only host/path plumbing and the sample loop were changed
  (SPEC.md §9). Behaviour-preservation quirks (SPEC.md §12) are retained on purpose.
- **`modules/local/*.nf`** — one process per computational unit, each with a `stub:` block.
- **`subworkflows/local/gas.nf`** — wires the GAS `ALL` routine (SPEC.md §3).
- **`workflows/wade.nf`** — parses the samplesheet and branches by organism.

> **Compatibility note:** `pairwiseAlignment()` moved from Biostrings to the `pwalign` package in
> Biostrings ≥ 2.77 (defunct in Biostrings). `common.R` calls it via a shim that prefers `pwalign`
> with a Biostrings fallback — algorithm and arguments are unchanged, so results are identical.

## Testing (SPEC.md §14)

```bash
nextflow run . -profile test,local_env -stub-run --outdir results_stub   # Layer 1: plumbing
nextflow run . -profile test,local_env           --outdir results        # Layer 2: env smoke
nf-test test --profile local_env                                         # Layer 3: regression snapshots

# Layer 4: scientific parity vs the original WADE R
Rscript tests/bin/gen_golden.R <wade_src> wade-data assets/test_data tests/fixtures/expected
python  tests/bin/compare_csv.py results/GAS/LabWareUpload_GAS_AMR.csv \
                                 tests/fixtures/expected/LabWareUpload_GAS_AMR.csv
```

`tests/bin/compare_csv.py` compares CSVs as *data* (sorts rows/cols, trims whitespace, normalises
NA) because WADE writes unquoted CSVs whose order can vary harmlessly.

`tests/fixtures/expected/` holds golden CSVs generated from the unmodified upstream R
(`tests/bin/gen_golden.R`) — the parity reference used by CI.

### Validation status (GAS)

- **Layers 1–4: green** (stub, end-to-end, nf-test snapshots, parity vs original R).
- **Layer 5 — determinism:** two independent fresh runs produce **identical** outputs; `-resume`
  re-caches all 9 tasks (0 resubmitted). Per-analysis FASTAs are collated with a content-based
  sort so re-runs are byte-identical (a `sort: true` collation ordered by random temp path and was
  fixed to `sort: { it.text }`).
- **Layer 5 — edge inputs:** a no-marker contig → clean `Susceptible` (all NEG); a present-but-empty
  contig completes without crashing the run (it yields spurious `POS/???` calls — the upstream R's
  behaviour when `blastout.txt` is empty; `Sample_Err` is reserved for a contig that can't be staged
  at all). A missing VCF → 23S count `NA`.

### Linting

`nextflow config -profile test,local_env` resolves cleanly. `nf-core pipelines lint` is **not** wired
up: it expects a pipeline generated from the nf-core template and crashes during discovery on this
hand-written layout. The pipeline follows nf-core *conventions* (module/subworkflow structure, meta
maps, `ext.args`/`modules.config`, `stub:` blocks, `nextflow_schema.json`, profiles); full template
lint would require scaffolding the template files (`CITATIONS.md`, the nf-core util subworkflows,
etc.) — tracked as future work. `nextflow lint` (the engine linter) needs Nextflow ≥ 25 (this repo
pins 24.10).

### Containers

`containers/environment.yml` (the runtime conda recipe: R + Biostrings + **pwalign** + BLAST+) is
validated by `-profile test,conda` (Nextflow builds it per-process and runs the pipeline in it).
`containers/Dockerfile` mirrors that package set on a pinned Bioconductor base; building/publishing
the image (for `-profile docker`/`singularity`) needs a Docker daemon and is left to CI/release.

## Adding the remaining organisms

Each org follows the GAS template:
1. Port its LabWare formatter into `bin/lib/labware/<org>_<test>.R` and add a `labware` dispatch
   branch in `bin/wade.R` (the per-sample `MASTERBLASTR`/`MLST`/`RRNA23S` modules are organism-agnostic and already done).
2. Add `modules/local/labware_<org>_<test>.nf`.
3. Write `subworkflows/local/<org>.nf` mirroring its row in SPEC.md §3.
4. Wire it into `workflows/wade.nf` (`ch_branched.<org>`).
5. Extend `gen_golden.R` + add fixtures + an nf-test.

## License

Apache-2.0, following upstream `phac-nml/wade`.
