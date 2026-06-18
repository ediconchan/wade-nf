# ============================================================================
# EMM -- single-sample port of R/EMM.R : EMM_pipeline()
# ----------------------------------------------------------------------------
# emm typing for GAS via blastn vs emm_trimmed.fasta (CDC trimmed emm DB).
# Sample loop removed; matching/quality logic kept verbatim (SPEC.md section 4e).
# Emits one row; COLLATE merges across samples.
# ============================================================================

run_emm <- function(org, sample, variable, contig, refdir, outdir, tempdir) {
  Org_id <- org
  CurrSampleNo <- sample
  Variable <- if (is.na(variable) || variable == "") NA else variable

  directorylist <- build_dirs(outdir, tempdir)
  reflist <- build_ref(refdir, directorylist, Variable)
  reflist <- prepare_lkup(reflist, directorylist)

  emmType <- ""
  emmSubtype <- ""
  emmTypeRep <- ""
  BP_ids <- ""
  emmComment <- ""

  CurrLocusLen <- 180L

  bad_emm_file <- paste0(reflist$Ref_Dir, "bad_emm_types.csv")
  df.bad_emm <- as_tibble(read.csv(bad_emm_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))

  if(!file.copy(contig, reflist$DestFile, overwrite = T))
  {
    emmType <- "Sample_Err"
    emmSubtype <- "Sample_Err"
    emmTypeRep <- "Sample_Err"
    BP_ids <- 0L
    emmComment <- "Sample not found"
  }

  if(emmType != "Sample_Err")
  {
    Blast_Out_File <- blastquery(directorylist, reflist, "emm_trimmed", 10e-50)

    info = file.info(Blast_Out_File)
    if(info$size == 0)
    {
      emmType <- "NF"
      emmSubtype <- "NF"
      emmTypeRep <- "NF"
      BP_ids <- 0L
      emmComment <- "No emm gene"
    }else
    {
      df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
      names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
      df.blastout <- distinct(df.blastout, Allele, .keep_all = TRUE)

      df.blastout_bad <- inner_join(df.blastout, df.bad_emm, by = "Allele")
      df.blastout_bad_100 <- filter(df.blastout_bad, Ident == 100 & Align == CurrLocusLen)
      emmComment <- tolower(paste0(df.blastout_bad_100$Allele, collapse = "/"))

      df.blastout_2 <- anti_join(df.blastout, df.bad_emm, by = "Allele")
      dfSize_blastout_2 <- nrow(df.blastout_2)

      if(dfSize_blastout_2 == 0)
      {
        df.blastout_2 <- df.blastout_bad_100
      }

      df.blastout_2$emmType <- sub("\\.\\d+","", df.blastout_2$Allele)

      df.blastout100 <- filter(df.blastout_2, Ident == 100 & Align == CurrLocusLen)

      if(nrow(df.blastout100 >= 1))
      {
        emmdf <- df.blastout100
      } else
      {
        emmdf <- dplyr::slice(df.blastout_2,1)
      }

      emmType <- tolower(paste0(emmdf$emmType, collapse = "/"))
      emmSubtype <- tolower(paste0(emmdf$Allele, collapse = "/"))
      BP_ids <- emmdf$Align[1] - emmdf$Mismatches[1]

      if(BP_ids <= 149 | is.na(BP_ids))
      {
        emmType <- "NF"
        emmSubtype <- "NF"
        emmTypeRep <- "NF"
        emmComment <- tolower(df.blastout_2$Allele[1])
        if(is.na(BP_ids) & nrow(df.blastout_bad > 0))
        {
          emmComment <- "partial match to bad allele"
        }
      }else
      {
        emmTypeRep[BP_ids >= 150] <- "NT"
        emmTypeRep[BP_ids == 180] <- emmSubtype
      }
    }
  }

  Output.df <- tibble(CurrSampleNo, emmType, emmTypeRep, emmSubtype, BP_ids, emmComment)
  names(Output.df) <- c("Sample", "Type", "Subtype_Rep", "Subtype", "bp_id", "Comments")

  cat(CurrSampleNo, "\t", emmType, "\t", emmTypeRep, "\t", emmSubtype, "\t", BP_ids, "\t", emmComment, "\n", sep = "")

  write.csv(Output.df, paste0(directorylist$output_dir, "output_profile_emm.csv"), quote = FALSE, row.names = FALSE)
  return(Output.df)
}
