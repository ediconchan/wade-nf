# ============================================================================
# LABWARE_GAS_TOXINS -- port of R/LabwareUpload_GAS_TOXINS.R
# ----------------------------------------------------------------------------
# Consumes the collated output_profile_GAS_TOXINS.csv -> LabWareUpload_GAS_TOXINS.csv
# ============================================================================

run_labware_gas_toxins <- function(profile, outdir) {
  outdir <- sub("/?$", "/", outdir)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  Output.df <- as_tibble(read.csv(profile, header = TRUE, sep = ",", stringsAsFactors = FALSE))

  Output.df$SampleProfile[Output.df$SampleProfile == ""] <- "Error"
  LabWare.df <- tibble(Output.df$SampleNo,
                       Output.df$speA_result,
                       Output.df$speC_result,
                       Output.df$speG_result,
                       Output.df$speH_result,
                       Output.df$speI_result,
                       Output.df$speJ_result,
                       Output.df$speK_result,
                       Output.df$speL_result,
                       Output.df$speM_result,
                       Output.df$smeZ_result,
                       Output.df$ssa_result,
                       Output.df$sagA_result,
                       Output.df$SampleProfile)

  names(LabWare.df) <- c("SampleNo", "speA", "speC", "speG", "speH", "speI",
                         "speJ", "speK", "speL", "speM","smeZ", "ssa", "sagA",
                         "TOXIN Profile")

  write.csv(LabWare.df, paste0(outdir, "LabWareUpload_GAS_TOXINS.csv"),
            quote = FALSE, row.names = FALSE)
  return(LabWare.df)
}
