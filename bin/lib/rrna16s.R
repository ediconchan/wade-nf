# ============================================================================
# RRNA16S_FORMAT -- port of R/16SrRNA.R : rRNA16S_pipeline()
# ----------------------------------------------------------------------------
# Post-processes the (collated) MASTERBLASTR rRNA16S profile into the LabWare
# species-ID table. <97% identity -> NF (SPEC.md section 4h). Runs once on the
# aggregated profile CSV (after COLLATE), so the sample loop is intrinsic to the
# table and unchanged.
# ============================================================================

run_rrna16s_format <- function(org, profile, outdir) {
  Org_id <- org
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outdir <- sub("/?$", "/", outdir)

  Output.df <- as_tibble(read.csv(profile, header = TRUE, sep = ",", stringsAsFactors = FALSE))

  Labware.df <- tibble(SampleNo = Output.df$SampleNo,
                       Allele = Output.df$X16S_allele,
                       Identification = Output.df$X16S_comments,
                       Percent_Match = Output.df$X16S_PctWT)

  Labware.df$Percent <- as.numeric(str_extract(Labware.df$Percent_Match,  "(?<=\\().+?(?=\\%)"))
  Labware.df$Allele[Labware.df$Allele == ""| Labware.df$Percent < 97] <- "NF"
  Labware.df$Identification[Labware.df$Identification == ""| Labware.df$Percent < 97] <- "NF"

  Labware.df <- select(Labware.df, -Percent)

  write.csv(Labware.df, paste0(outdir, "LabWareUpload_", Org_id, "_rRNA16S.csv"),
            quote = FALSE,  row.names = FALSE)
  return(Labware.df)
}
