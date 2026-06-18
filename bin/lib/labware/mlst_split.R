# ============================================================================
# MLST good/bad split -- post-collate step (R/MLST.R export section)
# ----------------------------------------------------------------------------
# Given the collated MLST/NGSTAR/NGMAST profile, writes:
#   LabWareUpload_<Org>_<Test>.csv          (all rows)
#   LabWareUpload_<Org>_<Test>_good.csv     (ST assigned)
#   LabWareUpload_<Org>_<Test>_bad.csv      (ST not found)
# NGMAST uses the "ST-NA" sentinel; others use NA (SPEC.md section 4d / 7).
# ============================================================================

run_mlst_split <- function(profile, org, test, outdir) {
  outdir <- sub("/?$", "/", outdir)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  Output <- as_tibble(read.csv(profile, header = TRUE, sep = ",", stringsAsFactors = FALSE))

  if(test == "NGMAST")
  {
    Output_good <- dplyr::filter(Output, ST != "ST-NA")
    Output_bad  <- dplyr::filter(Output, ST == "ST-NA")
  } else
  {
    Output_good <- filter(Output, !is.na(ST))
    Output_bad  <- filter(Output, is.na(ST))
  }
  write.csv(Output_good, paste0(outdir, "LabWareUpload_", org, "_", test, "_good.csv"), quote = FALSE, row.names = FALSE)
  write.csv(Output_bad,  paste0(outdir, "LabWareUpload_", org, "_", test, "_bad.csv"),  quote = FALSE, row.names = FALSE)
  write.csv(Output,      paste0(outdir, "LabWareUpload_", org, "_", test, ".csv"),      quote = FALSE, row.names = FALSE)
  return(Output)
}
