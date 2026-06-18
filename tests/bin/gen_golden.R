#!/usr/bin/env Rscript
# ============================================================================
# Generate golden fixtures by running the ORIGINAL phac-nml/wade R functions
# (sourced verbatim from a checkout) on the bundled example contigs, using the
# pinned wade-data. Outputs go to tests/fixtures/expected/<ORG>/ and are used by
# the Layer-4 parity test (SPEC.md section 14).
#
# pwalign is attached so the upstream (unmodified) pairwiseAlignment() call
# dispatches to it on Biostrings >= 2.77 -- the same shim the port uses.
#
# Usage: gen_golden.R <wade_src> <wade_data> <contigs_dir> <out_fixture_dir> <ORG>
#   ORG in {GAS, GBS, PNEUMO, GONO}.  contigs_dir must hold <ORG>.fasta (+ <ORG>.vcf
#   for PNEUMO/GONO 23S).
# ============================================================================
suppressWarnings(suppressMessages({
  library(plyr); library(tidyverse); library(tidyselect)
  library(data.table); library(Biostrings)
  if (requireNamespace("tidysq", quietly = TRUE)) library(tidysq)
  if (requireNamespace("pwalign", quietly = TRUE)) library(pwalign)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: gen_golden.R <wade_src> <wade_data> <contigs> <out_fixtures> <ORG>")
wade_src <- normalizePath(args[1]); wade_data <- normalizePath(args[2])
contigs  <- normalizePath(args[3]); fixtures_base <- args[4]; ORG <- args[5]
fixtures <- file.path(fixtures_base, ORG)
dir.create(fixtures, showWarnings = FALSE, recursive = TRUE)

for (f in list.files(file.path(wade_src, "R"), pattern = "\\.R$", full.names = TRUE)) {
  try(source(f), silent = TRUE)
}

# wade-data ships allele lookup FASTAs without BLAST indexes; index in place.
for (alleles in Sys.glob(file.path(wade_data, ORG, "*", "allele_lkup_dna", "*.fasta"))) {
  if (!file.exists(paste0(alleles, ".ndb")) && !file.exists(paste0(alleles, ".nin"))) {
    system(paste0("makeblastdb -in ", shQuote(alleles), " -dbtype nucl"))
  }
}

work <- file.path(tempdir(), paste0("orig_", ORG)); dir.create(work, showWarnings = FALSE, recursive = TRUE)
work <- paste0(normalizePath(work), "/")
dir.create(paste0(work, "Output"), showWarnings = FALSE)
dir.create(paste0(work, "temp"), showWarnings = FALSE)
writeLines(c("SampleNo,Variable", paste0(ORG, ",")), paste0(work, "list.csv"))
vcfdir <- if (ORG %in% c("PNEUMO", "GONO")) paste0(contigs, "/") else NA
dl <- data.frame(OrgID = ORG, LocalDir = work, SystemDir = paste0(wade_data, "/"),
                 ContigsDir = paste0(contigs, "/"), VCFDir = vcfdir)
write.csv(dl, paste0(work, "DirectoryLocations.csv"), row.names = FALSE, quote = FALSE)

# Run the org's ALL routine (SPEC.md section 3)
if (ORG == "GAS") {
  MASTER_pipeline("GAS", "AMR", "list", "list", work); labware_gas_amr("GAS", work)
  MASTER_pipeline("GAS", "TOXINS", "list", "list", work); labware_gas_toxins("GAS", work)
  MASTER_pipeline("GAS", "rRNA16S", "list", "list", work); rRNA16S_pipeline("GAS", work)
  MLST_pipeline("GAS", "MLST", "list", "list", work)
  EMM_pipeline("GAS", "list", work)
} else if (ORG == "GBS") {
  MASTER_pipeline("GBS", "AMR", "list", "list", work); labware_gbs_amr("GBS", work)
  MASTER_pipeline("GBS", "rRNA16S", "list", "list", work); rRNA16S_pipeline("GBS", work)
  MLST_pipeline("GBS", "MLST", "list", "list", work)
  SEROTYPE_pipeline("GBS", "list", "list", "SERO", work)
} else if (ORG == "PNEUMO") {
  MASTER_pipeline("PNEUMO", "AMR", "list", "list", work)
  rRNA23S_pipeline("PNEUMO", "list", work)
  labware_pneumo_amr("PNEUMO", work)
  MASTER_pipeline("PNEUMO", "VIRULENCE", "list", "list", work); labware_pneumo_virulence("PNEUMO", work)
  MLST_pipeline("PNEUMO", "MLST", "list", "list", work)
  SEROTYPE_pipeline("PNEUMO", "list", "list", "SERO", work)
} else if (ORG == "GONO") {
  MASTER_pipeline("GONO", "AMR_ALL", "list", "list", work)
  MLST_pipeline("GONO", "NGSTAR_AMR", "list", "list", work)
  rRNA23S_pipeline("GONO", "list", work)
  labware_gono_amr("GONO", work)
  MLST_pipeline("GONO", "MLST", "list", "list", work)
  MLST_pipeline("GONO", "NGSTAR", "list", "list", work)
  MLST_pipeline("GONO", "NGMAST", "list", "list", work)
}

for (n in list.files(paste0(work, "Output"), pattern = "\\.csv$")) {
  if (grepl("^LabWareUpload_|^output_profile_", n)) {
    file.copy(paste0(work, "Output/", n), file.path(fixtures, n), overwrite = TRUE)
  }
}
cat("\nGolden fixtures for ", ORG, " written to ", fixtures, "\n", sep = "")
