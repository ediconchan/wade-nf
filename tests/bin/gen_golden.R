#!/usr/bin/env Rscript
# ============================================================================
# Generate golden fixtures by running the ORIGINAL phac-nml/wade R functions
# (sourced verbatim from a checkout) on the bundled example contigs, using the
# pinned wade-data. Outputs are copied to tests/fixtures/expected/ and used by
# the Layer-4 parity test (SPEC.md section 14).
#
# pwalign is attached so the upstream (unmodified) pairwiseAlignment() call
# dispatches to it on Biostrings >= 2.77 -- the same compatibility shim the port
# uses, so the comparison stays apples-to-apples.
#
# Usage: gen_golden.R <wade_src_dir> <wade_data_dir> <contigs_dir> <out_fixture_dir>
# ============================================================================
suppressWarnings(suppressMessages({
  library(plyr); library(tidyverse); library(tidyselect)
  library(data.table); library(Biostrings)
  if (requireNamespace("pwalign", quietly = TRUE)) library(pwalign)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: gen_golden.R <wade_src> <wade_data> <contigs> <out_fixtures>")
wade_src <- normalizePath(args[1]); wade_data <- normalizePath(args[2])
contigs  <- normalizePath(args[3]); fixtures <- args[4]
dir.create(fixtures, showWarnings = FALSE, recursive = TRUE)

# source original logic
for (f in c("CommonFunctions.R", "MASTER.R", "MLST.R", "EMM.R", "16SrRNA.R",
            "LabwareUpload_GAS_AMR.R", "LabwareUpload_GAS_TOXINS.R")) {
  source(file.path(wade_src, "R", f))
}

# wade-data ships allele lookup FASTAs without BLAST indexes; the original
# pipeline expects them indexed (MakeblastDB_pipeline). Index in place so the
# upstream allele lookups behave as designed (deterministic; index-only).
for (alleles in Sys.glob(file.path(wade_data, "GAS", "*", "allele_lkup_dna", "*.fasta"))) {
  if (!file.exists(paste0(alleles, ".ndb")) && !file.exists(paste0(alleles, ".nin"))) {
    system(paste0("makeblastdb -in ", shQuote(alleles), " -dbtype nucl"))
  }
}

# build a WADE working dir + DirectoryLocations.csv (POSIX, trailing slashes)
work <- file.path(tempdir(), "orig_wade"); dir.create(work, showWarnings = FALSE, recursive = TRUE)
work <- paste0(normalizePath(work), "/")
dir.create(paste0(work, "Output"), showWarnings = FALSE)
dir.create(paste0(work, "temp"), showWarnings = FALSE)
writeLines(c("SampleNo,Variable", "GAS,"), paste0(work, "list.csv"))
dl <- data.frame(OrgID = "GAS", LocalDir = work,
                 SystemDir = paste0(wade_data, "/"),
                 ContigsDir = paste0(contigs, "/"), VCFDir = NA)
write.csv(dl, paste0(work, "DirectoryLocations.csv"), row.names = FALSE, quote = FALSE)

# run the original GAS routine
MASTER_pipeline("GAS", "AMR", "list", "list", work)
labware_gas_amr("GAS", work)
MASTER_pipeline("GAS", "TOXINS", "list", "list", work)
labware_gas_toxins("GAS", work)
MASTER_pipeline("GAS", "rRNA16S", "list", "list", work)
rRNA16S_pipeline("GAS", work)
MLST_pipeline("GAS", "MLST", "list", "list", work)
EMM_pipeline("GAS", "list", work)

# copy fixtures of interest
copy_if <- function(name) {
  src <- paste0(work, "Output/", name)
  if (file.exists(src)) file.copy(src, file.path(fixtures, name), overwrite = TRUE)
}
for (n in c("output_profile_GAS_AMR.csv", "LabWareUpload_GAS_AMR.csv",
            "output_profile_GAS_TOXINS.csv", "LabWareUpload_GAS_TOXINS.csv",
            "output_profile_GAS_rRNA16S.csv", "LabWareUpload_GAS_rRNA16S.csv",
            "output_profile_GAS_MLST.csv", "LabWareUpload_GAS_emm.csv")) copy_if(n)

cat("\nGolden fixtures written to ", fixtures, "\n", sep = "")
