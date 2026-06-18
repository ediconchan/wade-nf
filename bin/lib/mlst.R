# ============================================================================
# MLST -- single-sample port of R/MLST.R : MLST_pipeline()
# ----------------------------------------------------------------------------
# Sequence typing (MLST / NGSTAR / NGMAST) + NGSTAR_AMR mutation mode.
# Per-locus allele calls via locusblast()/mutationblast() (kept byte-for-byte in
# common.R, including the documented size-filter fall-through quirk -- SPEC.md
# section 12.1). Sample loop removed -> one row/sample; the ST left_join against
# profiles.txt is applied to that row (SPEC.md section 4d). NGMAST porB-/tbpB-/
# ST- prefixes and the good/bad split semantics are preserved.
# ============================================================================

run_mlst <- function(org, scheme, sample, variable, contig, refdir, outdir, tempdir) {
  Org_id <- org
  Test_id <- scheme
  SampleNo <- sample
  LocusID <- "list"
  Variable <- if (is.na(variable) || variable == "") NA else variable

  if(Test_id == "NGSTAR_AMR")
  {
    Test_id <- "NGSTAR"
    Test <- "NGSTAR_AMR"
  } else
  {
    Test <- "other"
  }

  directorylist <- build_dirs(outdir, tempdir)
  reflist <- build_ref(refdir, directorylist, Variable)
  reflist <- prepare_lkup(reflist, directorylist)
  unlink(paste0(directorylist$temp_dir, "blastout.txt"))

  blast_evalues.df <- as_tibble(read.csv(paste0(reflist$Ref_Dir, "blast_evalues.csv"),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Blast_evalue <- as.character(blast_evalues.df$contig[1])

  LocusList.df <- as_tibble(read.csv(reflist$Loci_List, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  locuslist <- LocusList.df$Locus_id
  NumLoci <- dim(LocusList.df)[1]

  # bad alleles (NGMAST only)
  df.bad_alleles <- NULL
  if(Test_id == "NGMAST")
  {
    df.bad_alleles <- as_tibble(read.csv(paste0(reflist$Ref_Dir, "bad_alleles.csv"),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE))
  }

  # profiles (not needed for NGSTAR_AMR mutation mode)
  if(Test == "other")
  {
    profiles.df <- as_tibble(read.csv(reflist$Profiles, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
    if(Test_id == "NGMAST")
    {
      names(profiles.df) <- c("ST", "porB", "tbpB")
    }
  }

  sampleinfo <- stage_contig(contig, SampleNo, Variable, reflist)

  if(Test == "NGSTAR_AMR")
  {
    allele_vec <- sapply(locuslist, mutationblast, sampleinfo = sampleinfo,
                         directorylist = directorylist, reflist = reflist,
                         Test_id = Test_id, df.bad_alleles = df.bad_alleles)
  } else
  {
    allele_vec <- sapply(LocusList.df$Locus_id, locusblast,
                         LocusList.df = LocusList.df, sampleinfo = sampleinfo,
                         reflist = reflist, Test_id = Test_id, Blast_evalue = Blast_evalue,
                         df.bad_alleles = df.bad_alleles, directorylist = directorylist)
  }

  # one-row allele table: SampleNo + per-locus allele
  Output.df <- tibble(SampleNo = SampleNo)
  for (loc in locuslist) Output.df[[loc]] <- as.character(allele_vec[[loc]])

  ##############################################################################
  # Lookup ST profiles
  ##############################################################################
  if(Test == "NGSTAR_AMR")
  {
    Output <- Output.df
  }else
  {
    if(Test_id == "NGSTAR")
    {
      profiles.df <- profiles.df %>% dplyr::rename("rRNA23S" = "X23S")
    }
    profiles.df <- profiles.df %>% dplyr::mutate(across(everything(), as.character))
    Output <- left_join(Output.df, profiles.df)
    if ("clonal_complex" %in% names(Output)) Output <- select(Output, -clonal_complex)
  }

  if(Test_id == "NGMAST")
  {
    Output$porB <- paste0("porB-", Output$porB)
    Output$tbpB <- paste0("tbpB-", Output$tbpB)
    Output$ST <- paste0("ST-", Output$ST)
  }

  if(Test == "NGSTAR_AMR")
  {
    outfile <- paste0(directorylist$output_dir, "output_profile_mut_", Org_id, "_NGSTAR.csv")
  } else
  {
    outfile <- paste0(directorylist$output_dir, "output_profile_", Org_id, "_", Test_id, ".csv")
  }
  write.csv(Output, outfile, quote = FALSE,  row.names = FALSE)
  return(Output)
}
