# ============================================================================
# LABWARE_PNEUMO_VIRULENCE -- port of R/LabwareUpload_PNEUMO_VIRULENCE.R
# ----------------------------------------------------------------------------
# Consumes the collated output_profile_PNEUMO_VIRULENCE.csv ->
# LabWareUpload_PNEUMO_VIRULENCE.csv (PI1/PI2/lytA/ply presence + rpoB).
# ============================================================================

run_labware_pneumo_virulence <- function(profile, outdir) {
  outdir <- sub("/?$", "/", outdir)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  Output.df <- as_tibble(read.csv(profile, header = TRUE, sep = ",", stringsAsFactors = FALSE))

  Output.df$SampleProfile[Output.df$SampleProfile == ""] <- "Err"

  LabWare.df <- tibble(Output.df$SampleNo,
                       Output.df$PI1_result,
                       Output.df$PI2_result,
                       Output.df$lytA_result,
                       Output.df$ply_result,
                       Output.df$rpoB_allele,
                       Output.df$rpoB_mutations,
                       Output.df$rpoB_PctWT,
                       Output.df$SampleProfile)
  names(LabWare.df) <- list("SampleNo", "PI1", "PI2", "lytA", "ply", "rpoB",
                            "rpoB_ID", "rpoB_match", "Profile")

  write.csv(LabWare.df, paste0(outdir, "LabWareUpload_PNEUMO_VIRULENCE.csv"),
            quote = FALSE, row.names = FALSE)
  return(LabWare.df)
}
