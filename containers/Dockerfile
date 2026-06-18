# ============================================================================
# wade-nf runtime image: R + Bioconductor + tidyverse + BLAST+ (SPEC.md section 10)
# ----------------------------------------------------------------------------
# Pinned Bioconductor release for reproducibility. Biostrings >= 2.77 moved
# pairwiseAlignment() to the pwalign package, which is installed here; the ported
# R uses it via a shim with a Biostrings fallback.
# ============================================================================
FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.source="https://github.com/phac-nml/wade-nf" \
      org.opencontainers.image.description="WADE molecular typing (Nextflow port): R + BLAST+" \
      org.opencontainers.image.licenses="Apache-2.0"

# Windows LMDB workaround; harmless on Linux (SPEC.md section 10)
ENV BLASTDB_LMDB_MAP_SIZE=1000000
ENV DEBIAN_FRONTEND=noninteractive

# NCBI BLAST+ (same source as upstream CI)
RUN apt-get update \
    && apt-get install -y --no-install-recommends ncbi-blast+ procps \
    && rm -rf /var/lib/apt/lists/*

# R packages: CRAN via install2.r (in the Bioconductor image) + Bioconductor
RUN install2.r --error --skipinstalled \
        plyr tidyverse tidyselect data.table optparse \
    && R -e 'BiocManager::install(c("Biostrings","pwalign"), update=FALSE, ask=FALSE)' \
    && R -e 'stopifnot(requireNamespace("pwalign", quietly=TRUE))'

CMD ["bash"]
