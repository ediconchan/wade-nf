# ============================================================================
# SEROTYPE -- single-sample port of R/SerotypeR.R : SEROTYPE_pipeline()
# ----------------------------------------------------------------------------
# 4-stage CPS serotyping for PNEUMO/GBS (SPEC.md section 4f):
#   1) serogroup via blastn vs reference_CPS.fasta
#   2) per-locus presence/absence (result), intact/pseudogene (pseudo),
#      serotype-determining AA substitutions (mutations), whole-gene allele match
#   3) progressively filter <serogroup>_loci_lookup.csv
#   4) remaining row(s) -> serotype
# Sample loop removed (one contig per task); locus loop preserved. Science kept
# verbatim. read_fasta()/tidysq replaced by base-R distinct (as in MASTERBLASTR).
# ============================================================================

run_serotype <- function(org, sample, contig, refdir, outdir, tempdir) {
  Org_id <- org
  SampleNo <- sample
  LocusID <- "list"
  Test_id <- "SERO"

  directorylist <- build_dirs(outdir, tempdir)
  reflist <- build_ref(refdir, directorylist, NA)
  reflist <- prepare_lkup(reflist, directorylist)

  CPS_types <- read.csv(paste0(reflist$Ref_Dir, "CPS_types_for_blast.csv"),
                        header = TRUE, sep = ",", stringsAsFactors = FALSE)
  CPS_types <- CPS_types$Serogroup

  blast_evalues.df <- as_tibble(read.csv(paste0(reflist$Ref_Dir, "blast_evalues.csv"),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Blast_evalue <- as.character(blast_evalues.df$contig[1])
  evalue_allele <- as.character(blast_evalues.df$allele[1])
  blast_wt_id_threshold <- as.numeric(blast_evalues.df$wt_id[1])

  sampleinfo <- stage_contig(contig, SampleNo, NA, reflist)

  SampleProfile <- ""
  serogroup <- ""
  serotype <- ""
  OutputLocus.df <- tibble(SampleNo = as.character(sampleinfo$CurrSampleNo))

  ##########################################################################
  # Get serogroup CPS locus for sample first
  ##########################################################################
  if (sampleinfo$Allele[1] == "Sample_Err")
  {
    serogroup <- "Sample_Err"
    serotype <- "Sample_Err"
    SampleProfile <- "Sample_Err"
  }else
  {
    BlastFormatCommand <- paste0("makeblastdb -in ", reflist$DestFile, " -dbtype nucl")
    try(system(BlastFormatCommand))

    Blast_Out_File <- blastquery(directorylist, reflist, "reference_CPS", Blast_evalue)
    info = file.info(Blast_Out_File)
    if(info$size == 0)
    {
      serotype <- "unknown"
      serogroup <- "unknown"
    }else
    {
      df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
      names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
      df.blastout <- arrange(df.blastout, dplyr::desc(bit))

      serotype1 <- df.blastout$Allele[1]
      SerotypeParts <- unlist(strsplit(serotype1, "_"))
      SerotypeParts <- SerotypeParts[1]
      serotype <- ifelse(substr(SerotypeParts, 1, 1) == "0", sub("^.", "", SerotypeParts), SerotypeParts)
      serotype[serotype %in% c("15B", "15C")] <- "15BC"

      if (Org_id == "PNEUMO")
      {
        serogroup <- as.character(substr(SerotypeParts, 1, 2))
        if(serogroup == "40") {serogroup <- "07"}
        if(serogroup == "44" || serogroup == "46"){serogroup <- "12"}
        if(serogroup == "38") {serogroup <- "25"}
        if(serogroup == "37") {serogroup <- "33"}
        if(serogroup == "42") {serogroup <- "35"}
      }
      if (Org_id == "GBS")
      {
        serogroup <- "GBS"
      }
      outfile_sero <- paste0(directorylist$output_dir, "output_blastout_serogroup.csv")
      write.csv(df.blastout, outfile_sero, quote = FALSE, row.names = F )
    }
  }

  cat("CPS Blast Serogroup ", serogroup, " - Serotype ", serotype, " !\n")

  if (serogroup == "unknown" | serogroup == "Sample_Err")
  {
    serotype <- paste("NT[", serogroup, "]", sep = "")
    SampleProfile <- serogroup
  }else
  {
    if(serotype %in% CPS_types)
    {
      SampleProfile <- "CPS operon type match"
      locus <- ""
      head(df.blastout, n = 5L)
    }else
    {
      LocusList <- paste0(reflist$Ref_Dir, serogroup, "_loci.csv")
      if(file.exists(LocusList))
      {
        LocusList.df <- read.csv(LocusList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
        SerotypeLookup.df <- read.csv(paste0(reflist$Ref_Dir, serogroup, "_loci_lookup.csv"),
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE,
                                      check.names = FALSE)
        NumLoci <- dim(LocusList.df)[1]
      }else
      {
        NumLoci <- 0
      }

      LocusMutationsFile <- paste0(reflist$Ref_Dir, serogroup, "_loci_mutations.csv")
      if(file.exists(LocusMutationsFile))
      {
        Mutns <- "Yes"
        LocusMutations.df <- read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
      } else
      {
        Mutns <- "No"
      }

      OutputLocus.df <- tibble(SampleNo = as.character(sampleinfo$CurrSampleNo))
      OutputLocusProfile.df <- OutputLocus.df

      p <- 1L
      for(p in 1L:NumLoci) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Locus Loop
      {
        ProfileSeparator <- ifelse(SampleProfile == "", "", ";")
        CurrLocus <- as.character(LocusList.df$Locus_id[p])
        locus <- gsub("GBS_", "", CurrLocus)
        locus_result_type <- as.character(LocusList.df$Result_type[p])

        LocusLkupDNA <- paste0(reflist$Lkup_Dir, locus, ".fasta")
        LocusLkupDNApresent <- ifelse(file.exists(LocusLkupDNA), TRUE, FALSE)

        resultcol <- paste0(locus, "_result")
        allelecol <- paste0(locus, "_allele")
        pseudocol <- paste0(locus, "_pseudo")
        mutationscol <- paste0(locus, "_mutations")
        OutputLocus.df[1, resultcol] <- NA
        OutputLocus.df[1, allelecol] <- NA
        OutputLocus.df[1, mutationscol] <- NA
        OutputLocus.df[1, pseudocol] <- NA

        LocusFile <- paste0(reflist$WT_Dir, locus, ".fasta")
        if(!file.exists(LocusFile))
        {
          sample_error.df <- tibble(Locus_ID = LocusID, Output = "Reference wildtype file not found.")
          OutputLocus.df <- OutputLocus.df %>% replace(is.na(.), "Locus_Err")
          return(sample_error.df)
        }else
        {
          BlastResult <- blasthits(reflist, directorylist, LocusFile, Blast_evalue)
          OutputLocus.df[1,2] <- BlastResult
          blastoutput <- readLines(paste0(directorylist$temp_dir, "blastout.txt"))

          if(locus_result_type == "result")
          {
            SampleProfile <- paste0(SampleProfile, ProfileSeparator, locus, "[", BlastResult, "]")
            cat(paste0(locus, "_result"), BlastResult, "\n", sep = "\t\t")

            if(LocusID == "list")
            {
              SerotypeLookup.df <- SerotypeLookup.df %>% filter(.[[resultcol]] == BlastResult|
                                                                 is.na(.[[resultcol]]))
            }
          }
          if((locus_result_type == "pseudo") & (BlastResult == "NEG"))
          {
            SampleProfile <- paste(SampleProfile, ProfileSeparator, locus, "[no pseudogene]", sep="")
            cat(pseudocol, "no gene present", "\n", sep = "\t\t")
          }
          if((locus_result_type == "mutations") & (BlastResult == "NEG"))
          {
            SampleProfile <- paste(SampleProfile, ProfileSeparator, locus, "[no mutn gene]", sep="")
            cat(mutationscol, "no gene present", "\n", sep = "\t\t")
          }
          if((locus_result_type == "allele") & (BlastResult == "NEG"))
          {
            SampleProfile <- paste(SampleProfile, ProfileSeparator, locus, "[no allele gene]", sep="")
            cat(allelecol, "no gene present", "\n", sep = "\t\t")
          }

          DNASeqLine_str <- ""
          WTDNASeqLine_str <- ""
          blastdetails <- ""
          if(BlastResult == "POS") #========================================== Positive BLAST Result Loop
          {
            blast_parsed <- parseblast(blastoutput, LocusID)

            WTDNASeqLine <- DNAString(blast_parsed$WTDNASeqLine_str)
            WTDNASeqLine_NoDash <- DNAString(blast_parsed$WTDNASeqLine_NoDash_str)
            DNASeqLine <- DNAString(blast_parsed$DNASeqLine_str)
            DNASeqLine_NoDash <- DNAString(blast_parsed$DNASeqLine_NoDash_str)

            blastAA <- AAconvert(WTDNASeqLine_NoDash, DNASeqLine_NoDash)

            if(locus_result_type == "pseudo")
            {
              WTend <- head(unlist(gregexpr('[*]', blastAA$WTAASeqLine_str[1])), n = 1)
              WTend <- ifelse(is.na(WTend), "0", WTend)
              queryend <- head(unlist(gregexpr('[*]', blastAA$AASeqLine_str[1])), n = 1)
              queryend <- ifelse(is.na(queryend), 0, queryend)
              diff <- ifelse(queryend == -1, 0, (WTend - queryend))
              lengthdiff <- blast_parsed$WTlength - blast_parsed$IDlength

              nonsense_mutation <- ifelse(diff > 6, "Disrupted",
                                          ifelse((blast_parsed$WTlength - blast_parsed$IDlength > 25),
                                          "Disrupted", "Intact"))

              OutputLocus.df[pseudocol] <- nonsense_mutation
              OutputLocusProfile.df[pseudocol] <- nonsense_mutation
              SampleProfile <- paste0(SampleProfile, ProfileSeparator, locus, "[", nonsense_mutation, "]")
              cat(pseudocol, nonsense_mutation, "\n", sep = "\t\t")
              if(LocusID == "list")
              {
                SerotypeLookup.df <- SerotypeLookup.df %>% filter(.[[pseudocol]] == nonsense_mutation|
                                                                  is.na(.[[pseudocol]]))
              }
            }

            motifs <- ""
            mutationspresent <- ""

            if(Mutns == "Yes")
            {
              motifs <- getmotifs(LocusMutations.df, Test_id, locus, LocusID, blastAA$AASeqLine_str)
            }

            OutputLocus.df[mutationscol] <- motifs
            if(locus_result_type == "mutations")
            {
              OutputLocusProfile.df[mutationscol] <- motifs
              SampleProfile <- paste(SampleProfile, ProfileSeparator, locus, "[", motifs, "]", sep="")
              cat(mutationscol, motifs, "\n", sep = "\t\t")
              if(LocusID == "list")
              {
                SerotypeLookup.df <- SerotypeLookup.df %>% filter(.[[mutationscol]] == motifs|
                                                                  is.na(.[[mutationscol]]))
              }
            }

            sink(reflist$DestFile, split=FALSE, append = FALSE)
            cat(">", sampleinfo$CurrSampleNo, "_", locus , "\n", blast_parsed$DNASeqLine_NoDash_str, sep ="")
            sink()

            if(locus_result_type == "allele")
            {
              Allele <- locusblast(locus, LocusList.df, sampleinfo, reflist, Test_id, Blast_evalue,
                                   directorylist = directorylist)
              OutputLocus.df[allelecol] <- Allele
              OutputLocusProfile.df[allelecol] <- Allele
              SampleProfile <- paste0(SampleProfile, ProfileSeparator, locus, "[", Allele, "]")
              cat(allelecol, Allele,"\n", sep = "\t\t")
              if(LocusID == "list")
              {
                SerotypeLookup.df <- SerotypeLookup.df %>% filter(.[[allelecol]] == Allele|
                                                                  is.na(.[[allelecol]]))
              }
            }
          }#================================================================== End BLAST positive
        }#******************************************************************** End WT gene found
      }#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End Locus Loop

      if(LocusID == "list")
      {
        NumResults <- dim(SerotypeLookup.df)[1]
        if(NumResults == 0)
        {
          serotype <- paste0("NT[", serotype, "]")
        }else
        {
          for(w in 1L:NumResults)
          {
            if(w == 1)
            {
              serotype <-  SerotypeLookup.df$Serotype[w]
            }else
            {
              serotype <- paste0(serotype, "/", SerotypeLookup.df$Serotype[w])
            }
          }
        }
      }
    }
  }

  SampleOutputProfile.df <- tibble(SampleNo = sampleinfo$CurrSampleNo, Serotype = serotype, Profile = SampleProfile)

  cat("\nProfile:  ", SampleProfile, sep = "")
  cat("\nSerotype: ", serotype, "\n", sep = "")

  outfile <- paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_SEROTYPE.csv")
  write.csv(SampleOutputProfile.df, outfile, quote = FALSE,  row.names = F)
  return(SampleOutputProfile.df)
}
