# ============================================================================
# LABWARE_GBS_AMR -- port of R/LabwareUpload_GBS_AMR.R : labware_gbs_amr()
# ----------------------------------------------------------------------------
# Consumes the collated output_profile_GBS_AMR.csv -> LabWareUpload_GBS_AMR.csv
# (+ _bad). MIC/interpretation logic kept verbatim. Note the preserved quirk
# (SPEC.md section 12): the bad-file filter references `amr_profile` (lowercase),
# which is NOT a column (the column is `AMR_profile`); dplyr data-masking falls
# back to the in-scope variable `amr_profile` (NA), so the bad file ends up empty.
# ============================================================================

run_labware_gbs_amr <- function(profile, outdir) {
  outdir <- sub("/?$", "/", outdir)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  Output.df <- as_tibble(read.csv(profile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  NumSamples <- dim(Output.df)[1]

  m <- 1
  for (m in 1L:NumSamples)
  {
    molec_profile <- NA
    allele_profile <- NA
    amr_profile <- NA
    lw_CurrSampleNo <- as.character(Output.df[m, "SampleNo"])

    if(Output.df[m,"SampleProfile"] == "Sample_Err")
    {
      sample_data.df <- tibble(lw_CurrSampleNo, lw_ermA = "Sample_Err",
                               lw_ermB = "Sample_Err", lw_ermT = "Sample_Err",
                               lw_mefAE = "Sample_Err", lw_lsaC = "Sample_Err",
                               lw_lnuB = "Sample_Err", lw_gyrA = "Sample_Err",
                               lw_parC = "Sample_Err", lw_tetM = "Sample_Err",
                               lw_tetO = "Sample_Err", lw_tetT = "Sample_Err",
                               lw_tetL = "Sample_Err", lw_cat = "Sample_Err",
                               lw_catQ = "Sample_Err", lw_pbp2x = "Sample_Err",
                               molec_profile = "Sample_Err",
                               ery_MIC = "Sample_Err", ery = "Sample_Err",
                               chl_MIC = "Sample_Err", chl = "Sample_Err",
                               lev_MIC = "Sample_Err", lev = "Sample_Err",
                               cli_MIC = "Sample_Err", cli = "Sample_Err",
                               tet_MIC = "Sample_Err", tet = "Sample_Err",
                               pen_MIC = "Sample_Err", pen = "Sample_Err",
                               AMR_profile = "Sample_Err", allele_profile = "Sample_Err")

    } else
    {
      ##### ermA #####
      lw_ermA <- as.character(Output.df[m, "ermA_result"])
      lw_ermAp_allele <- NA
      molec_profile_ermA <- NA

      if(lw_ermA == "POS")
      {
        lw_ermAp_result <- as.character(Output.df[m, "ermAp_result"])
        lw_ermAp_allele <- as.character(Output.df[m, "ermAp_allele"])
        lw_ermAp_mutation <- as.character(Output.df[m, "ermAp_mutations"])
        lw_ermAp_comment <- as.character(Output.df[m, "ermAp_comments"])
        lw_ermAp_motif <- as.character(Output.df[m, "ermAp_motifs"])

        if(lw_ermAp_motif == "WT/WT") {lw_ermAp_motif <- ""}
        if(lw_ermAp_result == "NEG")
        {
          molec_profile_ermA <- "ermA R"
          lw_ermAp_allele <- "ermAp 0"
        }else
        {
          if(lw_ermAp_motif == "D61G/P74Q")
          {
            molec_profile_ermA <- "ermA S D61G/P74Q"
            lw_ermAp_allele <- paste0("ermAp ", lw_ermAp_allele)
          }else
          {
            molec_profile_ermA <- paste("ermA R", lw_ermAp_motif, sep = " ")
            lw_ermAp_allele <-"ermAp 0"
          }
        }
      }

      ermB <- posneg_gene(Output.df, "ermB", m)
      ermT <- posneg_gene(Output.df, "ermT", m)
      mefAE <- posneg_gene(Output.df, "mefAE", m)
      lnuB <- posneg_gene(Output.df, "lnuB", m)
      lsaC <- posneg_gene(Output.df, "lsaC", m)
      gyrA <- allele1SNP(Output.df, "gyrA", "motifs", m)
      parC <- allele3SNPs(Output.df, "parC", "parC", "motifs", "D78", "S79", "D83", m)
      tetM <- posneg_gene(Output.df, "tetM", m)
      tetM$molec_profile_tetM[Output.df$tetM_comments[m] == "Disrupted"] <- "tetM Disrupted"
      tetO <- posneg_gene(Output.df, "tetO", m)
      tetT <- posneg_gene(Output.df, "tetT", m)
      tetL <- posneg_gene(Output.df, "tetL", m)
      cat <- posneg_gene(Output.df, "cat", m)
      catQ <- posneg_gene(Output.df, "catQ", m)
      dfrF <- posneg_gene(Output.df, "dfrF", m)
      dfrG <- posneg_gene(Output.df, "dfrG", m)
      pbp2x <- allele2SNPs(Output.df, "pbp2x", "pbp2x", "motifs", "400", "552", m)
      vanA <- posneg_gene(Output.df, "vanA", m)
      vanB <- posneg_gene(Output.df, "vanB", m)
      vanC <- posneg_gene(Output.df, "vanC", m)
      vanD <- posneg_gene(Output.df, "vanD", m)
      vanE <- posneg_gene(Output.df, "vanE", m)
      vanG <- posneg_gene(Output.df, "vanG", m)

      ######################## Now put them all together #######################
      molec_profile <- paste(molec_profile_ermA, ermB$molec_profile_ermB,
                             ermT$molec_profile_ermT, mefAE$molec_profile_mefAE,
                             lnuB$molec_profile_lnuB, lsaC$molec_profile_lsaC,
                             gyrA$molec_profile_gyrA, parC$molec_profile_parC,
                             tetM$molec_profile_tetM, tetO$molec_profile_tetO,
                             tetT$molec_profile_tetT, tetL$molec_profile_tetL,
                             cat$molec_profile_cat, catQ$molec_profile_catQ,
                             dfrF$molec_profile_dfrF, dfrG$molec_profile_dfrG,
                             pbp2x$molec_profile_pbp2x, vanA$molec_profile_vanA,
                             vanB$molec_profile_vanB, vanC$molec_profile_vanC,
                             vanD$molec_profile_vanD, vanE$molec_profile_vanE,
                             vanG$molec_profile_vanG, sep = "; ")
      molec_profile <- gsub("NA; ", "", molec_profile)
      molec_profile <- sub("; NA", "", molec_profile)

      if(molec_profile == "NA") {molec_profile <- "Wild Type"}

      allele_profile <- paste(lw_ermAp_allele, gyrA$allele_gyrA, parC$allele_parC,
                              pbp2x$allele_pbp2x, sep = ": ")
      allele_profile <- gsub("NA: ", "", allele_profile)
      allele_profile <- sub(": NA", "", allele_profile)

      ##### Erythromycin (ERY) #####
      if(str_detect(molec_profile, paste(c("ermA", "ermB", "ermT", "mefAE"),collapse = '|')))
      {
        ery_MIC <- ">= 2 ug/ml"; ery <- "Resistant"; ery_AMR <- "ERY-R"
      }else
      {
        ery_MIC <- "<= 0.25 ug/ml"; ery <- "Susceptible"; ery_AMR <- NA
      }

      ##### Clindamycin (CLI) #####
      if(str_detect(molec_profile, paste(c("ermB", "lsaC", "lnuB"), collapse = '|')))
      {
        cli_MIC <- ">= 1 ug/ml"; cli <- "Resistant"; cli_AMR <- "CLI-R"
      }else if(str_detect(molec_profile, "ermT"))
      {
        cli_MIC <- "<= 0.12 ug/ml"; cli <- "Inducible"; cli_AMR <- "CLI-Ind"
      }else
      {
        cli_MIC <- "<= 0.12 ug/ml"; cli <-"Susceptible"; cli_AMR <- NA
      }

      if(str_detect(molec_profile, "ermA S"))
      {
        cli_MIC <- "<= 0.12 ug/ml"; cli <- "Inducible"; cli_AMR <- "CLI-Ind"
      }
      if(str_detect(molec_profile, "ermA R"))
      {
        cli_MIC <- ">= 1 ug/ml"; cli <- "Resistant"; cli_AMR <- "CLI-R"
      }

      ##### Chloramphenicol (CHL) #####
      if(str_detect(molec_profile, paste(c("cat", "catQ"),collapse = '|')))
      {
        chl_MIC <- ">= 16 ug/ml"; chl <- "Resistant"; chl_AMR <- "CHL-R"
      }else
      {
        chl_MIC <- "<= 4 ug/ml"; chl <- "Susceptible"; chl_AMR <- NA
      }

      ##### Levofloxacin (LEV) #####
      if(gyrA$lw_gyrA == "Err" | parC$lw_parC == "Err")
      {
        lev_MIC <- "Error"; lev <- "Error"; lev_AMR <- "LEV-Err"
      }else if(str_detect(molec_profile, "gyrA S81"))
      {
        lev_MIC <- ">= 8 ug/ml"; lev <- "Resistant"; lev_AMR <- "LEV-R"
      }else if((str_detect(molec_profile, paste(c("S79", "D83"), collapse = '|'))))
      {
        lev_MIC <- "4 ug/ml"; lev <- "Intermediate"; lev_AMR <- "LEV-I"
      }else if(str_detect(molec_profile, "parC"))
      {
        lev_MIC <- "2 ug/ml"; lev <- "Susceptible"; lev_AMR <- NA
      }else
      {
        lev_MIC <- "<= 0.5 ug/ml"; lev <- "Susceptible"; lev_AMR <- NA
      }

      ##### Tetracycline (TET) #####
      if(str_detect(molec_profile, paste(c("tetM", "tetO", "tetT", "tetL"), collapse = '|')))
      {
        tet_MIC <- ">= 8 ug/ml"; tet <- "Resistant"; tet_AMR <- "TET-R"
      }else
      {
        tet_MIC <- "<= 1 ug/ml"; tet <- "Susceptible"; tet_AMR <- NA
      }
      if(str_detect(molec_profile, "tetM Disrupted"))
      {
        tet_MIC <- "<= 1 ug/ml"; tet <- "Susceptible"; tet_AMR <- NA
      }

      ##### Trimethoprim/Sulfamethoxazole (SXT) #####
      if(str_detect(molec_profile, paste(c("dfrG", "dfrF"), collapse = '|')))
      {
        sxt_MIC <- ">= 4/76 ug/ml"; sxt <- "Resistant"; sxt_AMR <- "SXT-R"
      }else
      {
        sxt_MIC <- "<= 0.5/9.5 ug/ml"; sxt <- "Susceptible"; sxt_AMR <- NA
      }

      ##### Penicillin (PEN) #####
      if(str_detect(pbp2x$lw_pbp2x, "Err"))
      {
        pen_MIC <- "Error"; pen <- "Error"; pen_AMR <- "PEN-Err"
      }else
      {
        if(str_detect(molec_profile, "pbp2x"))
        {
          pen_MIC <- ">= 0.12 ug/ml"; pen <- "Unknown"; pen_AMR <- "pen-DS"
        }else
        {
          pen_MIC <- "<= 0.12 ug/ml"; pen <- "Susceptible"; pen_AMR <- NA
        }
      }

      ##### Vancomycin (VAN) #####
      if(str_detect(molec_profile, paste(c("vanA", "vanB", "vanC", "vanD", "vanE", "vanG"), collapse = '|')))
      {
        van_MIC <- ">= 2 ug/ml"; van <- "Resistant"; van_AMR <- "VAN-R"
      }else
      {
        van_MIC <- "<= 1 ug/ml"; van <- "Susceptible"; van_AMR <- NA
      }

      ###################### Combine them for AMR profile ######################
      AMR_profile <- paste(ery_AMR, cli_AMR, chl_AMR, lev_AMR, tet_AMR, sxt_AMR,
                           pen_AMR, van_AMR, sep = "/")
      AMR_profile <- gsub("NA/", "", AMR_profile)
      AMR_profile <- sub("/NA", "", AMR_profile)

      if(AMR_profile == "NA") {AMR_profile <- "Susceptible"}

      sample_data.df <- tibble(lw_CurrSampleNo, lw_ermA, ermB$lw_ermB,
                               ermT$lw_ermT, mefAE$lw_mefAE, lsaC$lw_lsaC,
                               lnuB$lw_lnuB, gyrA$lw_gyrA, parC$lw_parC,
                               tetM$lw_tetM, tetO$lw_tetO, tetT$lw_tetT,
                               tetL$lw_tetL, cat$lw_cat, catQ$lw_catQ,
                               pbp2x$lw_pbp2x, molec_profile, ery_MIC, ery,
                               chl_MIC, chl, lev_MIC, lev, cli_MIC, cli,
                               tet_MIC, tet, pen_MIC, pen, AMR_profile,
                               allele_profile)
      sample_data.df <- sample_data.df %>% rename_at(vars(contains("$")), ~sub(".*\\$","",.))
    }

    if(m==1)
    {
      lw_Output.df <- tibble(sample_data.df)
    }else
    {
      lw_Output.df <- rbind(lw_Output.df, sample_data.df)
    }
  }

  # Preserved quirk: `amr_profile` (lowercase) is the in-scope NA variable, not a
  # column -> bad file ends up empty (see header note).
  lw_output_bad.df <- filter(lw_Output.df, amr_profile == "Error")

  write.csv(lw_Output.df, paste0(outdir, "LabWareUpload_GBS_AMR.csv"), quote = FALSE,  row.names = FALSE)
  write.csv(lw_output_bad.df, paste0(outdir, "LabWareUpload_GBS_AMR_bad.csv"), quote = FALSE,  row.names = FALSE)
  return(lw_Output.df)
}
