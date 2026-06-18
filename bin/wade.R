#!/usr/bin/env Rscript
# ============================================================================
# wade.R -- single-sample, path-explicit CLI for the WADE analysis package.
# ----------------------------------------------------------------------------
# Dispatches one analysis on one sample with explicit paths, replacing the
# Shiny app / DirectoryLocations.csv / list.csv host coupling (SPEC.md section 9).
# Each subcommand sources the refactored WADE functions from ./lib and runs them.
#
#   wade.R masterblastr --org GAS --test AMR  --sample S1 --contig in.fasta --refdir REF/ --outdir .
#   wade.R mlst         --org GONO --scheme NGSTAR --sample S1 --contig in.fasta --refdir REF/ --outdir .
#   wade.R emm          --org GAS --sample S1 --contig in.fasta --refdir REF/ --outdir .
#   wade.R rrna23s      --org PNEUMO --sample S1 --vcf in.vcf --outdir .
#   wade.R rrna16s      --org GAS --profile output_profile_GAS_rRNA16S.csv --outdir .
#   wade.R labware      --org GAS --test AMR --profile output_profile_GAS_AMR.csv --outdir .
#   wade.R mlst_split   --org GONO --test NGMAST --profile collated.csv --outdir .
# ============================================================================

suppressWarnings(suppressMessages(library(optparse)))

# --- locate this script's directory so we can source ./lib/* -----------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  env <- Sys.getenv("WADE_LIB", unset = NA)
  if (!is.na(env)) return(env)
  return(getwd())
}
SCRIPT_DIR <- get_script_dir()
LIB <- file.path(SCRIPT_DIR, "lib")

source(file.path(LIB, "common.R"))

# --- subcommand + remaining args ---------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 1) {
  stop("Usage: wade.R <masterblastr|mlst|emm|rrna23s|rrna16s|labware|mlst_split> [options]")
}
subcmd <- argv[1]
rest   <- argv[-1]

opt_list <- list(
  make_option("--org",      type = "character", default = NA),
  make_option("--test",     type = "character", default = NA),
  make_option("--scheme",   type = "character", default = NA),
  make_option("--sample",   type = "character", default = NA),
  make_option("--variable", type = "character", default = NA),
  make_option("--contig",   type = "character", default = NA),
  make_option("--vcf",      type = "character", default = NA),
  make_option("--refdir",   type = "character", default = NA),
  make_option("--profile",  type = "character", default = NA),
  make_option("--outdir",   type = "character", default = "."),
  make_option("--tempdir",  type = "character", default = NA)
)
opt <- parse_args(OptionParser(option_list = opt_list), args = rest)

tempdir_or_default <- function(o) if (is.na(o$tempdir)) tempfile("wade_") else o$tempdir

switch(subcmd,
  "masterblastr" = {
    source(file.path(LIB, "masterblastr.R"))
    run_masterblastr(opt$org, opt$test, opt$sample, opt$variable, opt$contig,
                     opt$refdir, opt$outdir, tempdir_or_default(opt))
  },
  "mlst" = {
    source(file.path(LIB, "mlst.R"))
    run_mlst(opt$org, opt$scheme, opt$sample, opt$variable, opt$contig,
             opt$refdir, opt$outdir, tempdir_or_default(opt))
  },
  "emm" = {
    source(file.path(LIB, "emm.R"))
    run_emm(opt$org, opt$sample, opt$variable, opt$contig,
            opt$refdir, opt$outdir, tempdir_or_default(opt))
  },
  "rrna23s" = {
    source(file.path(LIB, "rrna23s.R"))
    run_rrna23s(opt$org, opt$sample, opt$vcf, opt$outdir)
  },
  "rrna16s" = {
    source(file.path(LIB, "rrna16s.R"))
    run_rrna16s_format(opt$org, opt$profile, opt$outdir)
  },
  "labware" = {
    test <- opt$test
    if (opt$org == "GAS" && test == "AMR") {
      source(file.path(LIB, "labware/gas_amr.R"));     run_labware_gas_amr(opt$profile, opt$outdir)
    } else if (opt$org == "GAS" && test == "TOXINS") {
      source(file.path(LIB, "labware/gas_toxins.R"));  run_labware_gas_toxins(opt$profile, opt$outdir)
    } else {
      stop(paste0("No labware formatter for org=", opt$org, " test=", test))
    }
  },
  "mlst_split" = {
    source(file.path(LIB, "labware/mlst_split.R"))
    run_mlst_split(opt$profile, opt$org, opt$test, opt$outdir)
  },
  stop(paste0("Unknown subcommand: ", subcmd))
)
