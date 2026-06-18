# ============================================================================
# RRNA23S -- single-sample port of R/23SrRNA.R : rRNA23S_pipeline()
# ----------------------------------------------------------------------------
# Counts mutated 23S rRNA alleles (0-4) from a FreeBayes VCF. No BLAST.
# Organism-specific grep strings and the AO/DP -> allele-fraction bins are kept
# literally (SPEC.md section 12.7). The original per-line loop resets the count
# on every non-matching line, so the LAST line processed wins -- this quirk is
# preserved exactly. Missing VCF -> NA (SPEC.md section 4g).
# ============================================================================

run_rrna23s <- function(org, sample, vcf, outdir) {

  Org_id <- org
  CurrSampleNo <- sample

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outdir <- sub("/?$", "/", outdir)

  if (Org_id == "GONO")
  {
    rRNA23S_position1 <- "23S4_NCCP11945\t2045"
    rRNA23S_position2 <- "23S4_NCCP11945\t2597"
  }
  if (Org_id == "PNEUMO")
  {
    rRNA23S_position1 <- "23S_rRNA_R6_sprr02\t2061"
    rRNA23S_position2 <- "23S_rRNA_R6_sprr02\t2613"
  }

  A2059G <- 99L
  C2611T <- 99L

  QueryFile <- vcf
  if (is.null(vcf) || is.na(vcf) || vcf == "" || !file.exists(vcf))
  {
    A2059G <- NA
    C2611T <- NA
    QueryFile <- "File not found"
    cat("Sample Number ", CurrSampleNo, " not found!\n", sep = "")
  }

  if(QueryFile != "File not found")
  {
    allele_fraction <- 0L
    con <- file(QueryFile, open="r")
    linn <- readLines(con)
    close(con)
    for (i in 1:length(linn))
    {
      if(str_detect(linn[i], rRNA23S_position1))  #2597 & 2045 for GONO; 2061 for pneumo
      {
        vcf_parts <- strsplit(linn[i], ";")
        vcf_parts <- unlist(vcf_parts)
        DP1 <- vcf_parts[8]
        DP<-as.integer(substr(DP1, 4, 10))
        AO1 <- vcf_parts[6]
        AO<-as.integer(substr(AO1, 4, 10))

        allele_fraction <- as.integer((AO/DP)*100)
        if(allele_fraction < 14){A2059G <- 0L}
        if((allele_fraction >= 15) && (allele_fraction <= 34)) {A2059G <- 1L}
        if((allele_fraction >= 35) && (allele_fraction <= 64)) {A2059G <- 2L}
        if((allele_fraction >= 65) && (allele_fraction <= 84)) {A2059G <- 3L}
        if(allele_fraction >= 85) {A2059G <- 4L}
      }else
      {
        A2059G <- 0L
      }

      if(str_detect(linn[i], rRNA23S_position2))  #2597 & 2045 for GONO; 2061 for pneumo
      {
        vcf_parts2 <- strsplit(linn[i], ";")
        vcf_parts2 <- unlist(vcf_parts2)
        DP2_str <- vcf_parts2[8]
        DP2 <- as.integer(substr(DP2_str, 4, 10))   #total number of reads (depth of reads)
        AO2_str <- vcf_parts2[6]
        AO2 <- as.integer(substr(AO2_str, 4, 10))  #alternate observations from reference

        allele_fraction <- as.integer((AO2/DP2)*100)
        if(allele_fraction < 14){C2611T <- 0L}
        if((allele_fraction >= 15) && (allele_fraction <= 34)) {C2611T <- 1L}
        if((allele_fraction >= 35) && (allele_fraction <= 64)) {C2611T <- 2L}
        if((allele_fraction >= 65) && (allele_fraction <= 84)) {C2611T <- 3L}
        if(allele_fraction >= 85) {C2611T <- 4L}
      }else
      {
        C2611T <- 0L
      }
    }
  }

  cat(CurrSampleNo, "\tAC_2059: ", A2059G, "\tAC_2611: ", C2611T,"\n", sep = "")

  OutputProfile.df <- tibble(CurrSampleNo, A2059G, C2611T)
  names(OutputProfile.df) <- c("SampleNo", "A2059G", "C2611T")

  outfile <- paste0(outdir, "output_profile_", Org_id, "_rRNA23S.csv")
  write.csv(OutputProfile.df, outfile, row.names = F)

  return(OutputProfile.df)
}
