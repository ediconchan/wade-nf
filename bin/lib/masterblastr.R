# ============================================================================
# MASTERBLASTR -- single-sample port of R/MASTER.R : MASTER_pipeline()
# ----------------------------------------------------------------------------
# Presence/absence + allele + mutation calling per locus (SPEC.md section 4c).
# The original looped over a list.csv of samples and over loci.csv; here the
# sample loop is removed (Nextflow scatters one contig per task) and the locus
# loop is preserved inside (SPEC.md section 9.3). The per-locus science is kept
# byte-for-byte. Behaviour-preservation quirks (SPEC.md section 12): disruption
# rule (* or coverage<80), wt_id NEG threshold, porB1a/1b, BLOSUM50 gapOpen -2
# gapExtend -8 alignment, trailing-* trim -- all retained.
# ============================================================================

run_masterblastr <- function(org, test, sample, variable, contig, refdir, outdir, tempdir) {

  Org_id <- org
  Test_id <- test
  SampleNo <- sample
  LocusID <- "list"
  Variable <- if (is.na(variable) || variable == "") NA else variable

  # Test_id remaps preserved from MASTER_pipeline (SPEC.md section 3 note)
  if (Test_id == "AMR_ALL") Test_id <- "AMR"
  if (Org_id == "PNEUMO" && Test_id %in% c("AMR_ALL", "AMR_2")) Test_id <- "AMR"

  directorylist <- build_dirs(outdir, tempdir)
  reflist <- build_ref(refdir, directorylist, Variable)

  dna_file    <- paste0(directorylist$output_dir, "output_dna.fasta")
  dna_nf_file <- paste0(directorylist$output_dir, "output_dna_notfound.fasta")
  aa_file     <- paste0(directorylist$output_dir, "output_aa.fasta")
  # fresh per-task FASTAs (original cleared these via removefiles())
  unlink(c(dna_file, dna_nf_file, aa_file,
           paste0(directorylist$output_dir, "output_dna_notfound_distinct.fasta")))

  #-----------------------------------------------------------------------------
  # Blast Evalues
  blast_evalues.df <- as_tibble(read.csv(paste0(reflist$Ref_Dir, "blast_evalues.csv"),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Blast_evalue <- as.character(blast_evalues.df$contig[1])
  evalue_allele <- as.character(blast_evalues.df$allele[1])
  blast_wt_id_threshold <- as.numeric(blast_evalues.df$wt_id[1])
  #-----------------------------------------------------------------------------

  # Single sample (no list.csv)
  SampleList.df <- tibble(SampleNo, Variable)

  # Locus list from reference loci.csv
  LocusList.df <- as_tibble(read.csv(reflist$Loci_List, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  NumLoci <- dim(LocusList.df)[1]

  LocusMutationsFile <- paste0(reflist$Ref_Dir, "loci_mutations.csv")
  if(file.exists(LocusMutationsFile))
  {
    LocusMutations.df <- as_tibble(read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
    Mutns <- "Yes"
  } else
  {
    Mutns <- "No"
  }

  # Copy + index allele lookup DBs into a writable temp dir
  reflist <- prepare_lkup(reflist, directorylist)

  m <- 1L
  sampleinfo <- stage_contig(contig, SampleNo, Variable, reflist)
  CurrSample.df <- filter(SampleList.df, SampleNo == sampleinfo$CurrSampleNo)
  SampleProfile <- ""

  if(sampleinfo$Allele[1] != "Sample_Err"){
    #make blast database of contig file
    BlastFormatCommand <- paste0("makeblastdb -in ", reflist$DestFile, " -dbtype nucl")
    try(system(BlastFormatCommand))
  }

  p <- 1L
  for(p in 1L:NumLoci)   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start Locus Loop
  {
    AlleleLine <- ""
    AlleleInfo <- ""
    AlleleInfo[1:4] <- ""
    blast_parsed <- ""
    blastAA <- ""
    IDpercent2 <- ""
    BlastResult <- NA
    motifs <- ""

    CurrLocus <- as.character(LocusList.df[p,1])
    LocusLkupDNA <- paste0(reflist$Lkup_Dir, CurrLocus, ".fasta")
    LocusLkupDNApresent <- ifelse(file.exists(LocusLkupDNA), TRUE, FALSE)

    LocusFile <- paste0(reflist$WT_Dir, CurrLocus, ".fasta")
    if(!file.exists(LocusFile))
    {
      sample_error.df <- tibble(Locus_ID = CurrLocus, Output = "File Not Found", stringsAsFactors = FALSE)
      return(sample_error.df)
    }

    AlleleInfo[1][sampleinfo$Allele == "Sample_Err"] <- "Sample_Err"
    if(AlleleInfo[1] != "Sample_Err") #....................................... Not Sample Error
    {
      BlastResult <- blasthits(reflist, directorylist, LocusFile, Blast_evalue)
      AlleleInfo[1] <- ifelse(BlastResult == "POS", "POS", "NEG")

      if(BlastResult == "POS") #============================================== Positive BLAST Result Loop
      {
        blastoutput <- readLines(paste0(directorylist$temp_dir, "blastout.txt"))
        blast_parsed <- parseblast(blastoutput, LocusID)

        WTDNASeqLine <- DNAString(blast_parsed$WTDNASeqLine_str[1])
        WTDNASeqLine_NoDash <- DNAString(blast_parsed$WTDNASeqLine_NoDash_str[1])
        DNASeqLine <- DNAString(blast_parsed$DNASeqLine_str[1])
        DNASeqLine_NoDash <- DNAString(blast_parsed$DNASeqLine_NoDash_str[1])

        blastAA <- AAconvert(WTDNASeqLine_NoDash, DNASeqLine_NoDash)

        WTAASeqLine <- AAString(blastAA$WTAASeqLine_str)
        AASeqLine <- AAString(blastAA$AASeqLine_str)
        globalAlign_AA <- pairwise_align(AASeqLine, WTAASeqLine, substitutionMatrix = "BLOSUM50",
                                            gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
        WTAASeqLine_aln <- aln_subject(globalAlign_AA)
        WTAASeqLine_aln_str <- toString(WTAASeqLine_aln)
        AASeqLine_aln <- aln_pattern(globalAlign_AA)
        AASeqLine_aln_str <- toString(AASeqLine_aln)

        # DNA "alignment" sequence that only shows mutations
        DNASeqLine_aln <- ""
        for(j in 1:str_length(blast_parsed$WTDNASeqLine_str))
        {
          if(str_sub(blast_parsed$WTDNASeqLine_str, j, j) == str_sub(blast_parsed$DNASeqLine_str, j, j))
          {
            DNASeqLine_aln <- paste0(DNASeqLine_aln, ".")
          } else
          {
            DNASeqLine_aln <- paste0(DNASeqLine_aln, str_sub(blast_parsed$DNASeqLine_str, j, j))
          }
        }

        # Mutations
        mutationspresent <- ""
        AASeqLine_aln_disp <- ""
        mutations <- ""

        for(k in 1:str_length(WTAASeqLine_aln_str))
        {
          if(str_sub(WTAASeqLine_aln_str, k, k) == str_sub(AASeqLine_aln_str, k, k))
          {
            AASeqLine_aln_disp <- paste0(AASeqLine_aln_disp, ".")
          } else
          {
            AASeqLine_aln_disp <- paste0(AASeqLine_aln_disp, str_sub(blastAA$AASeqLine_str, k, k))
            mutations <- paste0(mutations, str_sub(WTAASeqLine_aln_str, k, k), k, str_sub(blastAA$AASeqLine_str, k, k), " ")
          }
        }

        if (Mutns == "Yes")
        {
          motifs <- getmotifs(LocusMutations.df, Test_id, CurrLocus, LocusID, AASeqLine_aln_str, AASeqLine_aln_disp)
        }

        #Check for disrupted genes
        AlleleInfo[4] <- ifelse(str_detect(AASeqLine_aln_disp, "[*]"), "Disrupted",
                                ifelse(blast_parsed$IDcoverage < 80, "Disrupted", ""))

        # Lookup DNA alleles
        Seq_File <- paste0(directorylist$temp_dir, "querygene.fasta")
        sink(Seq_File, split=FALSE, append = FALSE)
        cat(">", sampleinfo$CurrSampleNo, "_", CurrLocus , "\n", blast_parsed$DNASeqLine_NoDash_str, sep ="")
        sink()

        ExactMatchFound <- FALSE

        if(LocusLkupDNApresent) #********************************************* Locus Lookup DNA Loop
        {
          BlastCommand <- paste0("blastn -query ", Seq_File,
                                 " -db ", LocusLkupDNA,
                                 " -out ", directorylist$temp_dir, "blastout2.txt",
                                 " -evalue ", evalue_allele,
                                 " -num_alignments 1")
          system(BlastCommand, intern = TRUE)

          blastoutput2 <- readLines(paste0(directorylist$temp_dir, "blastout2.txt"))

          querylength <- grep("Length=", blastoutput2, value = TRUE)
          querylength <- sub("Length=", "", querylength)
          lengthmatch <- ifelse(querylength[1] == querylength[2], TRUE, FALSE)
          IDLine2 <- grep("Identities", blastoutput2, value = TRUE)
          ExactMatchFound <- sum(str_detect(IDLine2, "(100%)"))

          if(lengthmatch == TRUE & ExactMatchFound > 0)
          {
            AlleleLine <- grep(">", blastoutput2, value = TRUE)
            AlleleParts <- unlist(strsplit(AlleleLine, "_"))
            AlleleInfo[2] <- AlleleParts[2]
            AlleleInfo[3] <- AlleleParts[3]
            AlleleInfo[4] <- AlleleParts[4]
          } else
          {
            AlleleInfo[2] <- "NF"
            AlleleInfo[3] <- "???"
            AlleleInfo[4] <- ""
          }
        } #******************************************************************* End Locus Lookup DNA Loop

        # Create Sample Output
        sink(dna_file, split=FALSE, append = TRUE)
        cat(">", CurrLocus, "_", sampleinfo$CurrSampleNo, "_", CurrLocus, AlleleInfo[2], "_",
            AlleleInfo[3], "_", sampleinfo$CurrSampleVar, "_", motifs, "\n", blast_parsed$DNASeqLine_NoDash_str, "\n", sep ="")
        sink()

        if(AlleleInfo[2] == "NF")
        {
          sink(dna_nf_file, split=FALSE, append = TRUE)
          cat(">", CurrLocus, "_", motifs, "_", sampleinfo$CurrSampleNo, "_",
              sampleinfo$CurrSampleVar, "_", "\n", blast_parsed$DNASeqLine_NoDash_str, "\n", sep ="")
          sink()
        }

        sink(aa_file, split=FALSE, append = TRUE)
        cat(">", CurrLocus, "_", sampleinfo$CurrSampleNo, "_",
            sampleinfo$CurrSampleVar, "\n", blastAA$AASeqLine_str, "\n", sep ="")
        sink()

        IDpercent2 <- paste0(blast_parsed$IDlength, "/", blast_parsed$WTlength,
                            " (", blast_parsed$IDpercent, ")")

        # special notes for porB or poor quality BLAST
        AlleleInfo[4][CurrLocus == "porB" & blast_parsed$IDcoverage < 90] <- "porB1a"
        AlleleInfo[4][CurrLocus == "porB" & blast_parsed$IDcoverage >= 90] <- "porB1b"
        AlleleInfo[1][blast_parsed$IDcoverage <= blast_wt_id_threshold] <- "NEG"
      } else #================================================================ End Positive BLAST Result Loop
      {
        blast_parsed <- tibble(IDLine = "",
                               IDLine_trimmed = "")
        IDpercent2 <- ""
        motifs <- ""
      }
    } else #.................................................................. Close Not Sample Error Loop
    {
      AlleleInfo[2] <- "Sample_Err"
      AlleleInfo[3] <- "Sample_Err"
      AlleleInfo[4] <- "Sample_Err"
      blast_parsed <- tibble(IDLine = "Sample_Err",
                              IDLine_trimmed = "Sample_Err")
      IDpercent2 <- "Sample_Err"
      motifs <- "Sample_Err"
    }

    headers <- c("result", "allele", "mutations", "comments", "BlastID",
                 "PctWT", "motifs")
    headers <- paste(CurrLocus, headers, sep = "_")

    if(p==1)
    {
      OutputLocus.df <- tibble(AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], blast_parsed$IDLine_trimmed, IDpercent2, motifs)
      names(OutputLocus.df) <- headers
    }else
    {
      OutputLocus2.df <- tibble(AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], blast_parsed$IDLine_trimmed, IDpercent2, motifs)
      names(OutputLocus2.df) <- headers
      OutputLocus.df <- cbind(OutputLocus.df, OutputLocus2.df)
    }

    # Output Profile
    ProfileEntry <- ""
    if(AlleleInfo[1] == "POS")
    {
      ProfileEntry <- CurrLocus
      if(LocusLkupDNApresent == TRUE)
      {
        ProfileEntry <- paste(CurrLocus, AlleleInfo[3], sep = " ")
        ProfileEntry[AlleleInfo[3] %in% c("WT", "WT/WT", "WT/WT/WT", "WT/WT/WT/WT")] <- ""
      }
      SampleProfile[SampleProfile != "" & ProfileEntry != ""] <- paste(SampleProfile, ProfileEntry, sep = "-")
      SampleProfile[SampleProfile == ""] <- ProfileEntry
    }
    SampleProfile[AlleleInfo[1] == "Sample_Err"] <- "Sample_Err"

    cat(sampleinfo$CurrSampleNo, CurrLocus, AlleleInfo[1], AlleleInfo[2], AlleleInfo[3],
        AlleleInfo[4], blast_parsed$IDLine_trimmed, IDpercent2, motifs, "\n", sep = "\t")

  } #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End Locus List Loop

  SampleProfile <- ifelse(SampleProfile == "", "Susceptible", SampleProfile)
  SampleProfile.df <- tibble(SampleProfile)
  OutputLocus1.df <- cbind(CurrSample.df, OutputLocus.df, SampleProfile.df)

  OutputProfile.df <- tibble(OutputLocus1.df)
  names(OutputProfile.df) <- names(OutputLocus1.df)

  outfile <- paste0(directorylist$output_dir, "output_profile_", Org_id, "_", Test_id, ".csv")
  write.csv(OutputProfile.df, outfile, quote = FALSE,  row.names = F)

  # de-dup not-found alleles (original used tidysq::read_fasta; do it with base R
  # to avoid the extra dependency while preserving distinct-by-sequence output)
  if(file.exists(dna_nf_file))
  {
    nf_lines <- readLines(dna_nf_file)
    hdr_idx <- grep("^>", nf_lines)
    if (length(hdr_idx) > 0) {
      seqs <- nf_lines[hdr_idx + 1]
      keep <- !duplicated(seqs)
      out <- character(0)
      for (i in which(keep)) out <- c(out, nf_lines[hdr_idx[i]], seqs[i])
      unique_nf_file <- paste0(directorylist$output_dir, "output_dna_notfound_distinct.fasta")
      writeLines(out, unique_nf_file)
    }
  }

  cat("\n\nDone! output_profile_", Org_id, "_",  Test_id, ".csv written\n", sep = "")
  return(OutputProfile.df)
}
