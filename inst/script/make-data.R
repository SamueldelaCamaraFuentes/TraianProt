# DATA GENERATION OVERVIEW
# ---------------------------------------
# This script describes the generation/source of the datasets included 
# in the 'inst/extdata' directory of the TraianProt package.
#
# Author: Samuel de la Cámara
# Date: 2026-01-26


# Data documentation
# ---------------------------------------
# File: metadata_MaxQuant.tsv
# Source: Lab data format 
# License: Free for use for package demonstration (CC-BY / Public Domain)

# PROCEDURE:
# 1. Metadata file was generated following instructions for MaxQuant metadata included in Traianprot´s tutorial: "tutorial.pdf"

# Note: This file is static and is not generated via code at runtime.

# ---------------------------------------

# File: proteinGroups.txt.gz (compressed)
# Source: Internal laboratory data (Anonymized) / [PRIDE ID: PXD040804]
# License: CC-BY / Public Domain

# PROCEDURE:
# 1. Raw files from PXD040804 were re-processed using MaxQuant v2.1.4.
# 2. The Candida Genome Database (release 2020_06, 6209 sequences) was used.
# 3. Search parameters included Carbamidomethylation of cysteines as fixed modification, oxidation of methionine and N-terminal acetylation as variable modifications and Trypsin/P as proteolytic enzyme with a maximum of 2 missed cleavages allowed. The precursor mass tolerance was 10 ppm and the fragment mass tolerance was 0.02 Da.
# 4. The protein group expression table was ompressed using gzip.


# ==============================================================================
# DATASET 2: proteinGroups.txt
# ==============================================================================
# File: proteinGroups.txt.gz
# Source: Output from MaxQuant v2.x processing of the raw files associated 
#         with the metadata above.
# License: CC-BY / Public Domain

# PROCEDURE:
# 1. Raw mass spectrometry files were processed using MaxQuant with standard settings.
# 2. To reduce package size and computation time for vignettes, the original 
#    'proteinGroups.txt' file (which was >100MB) was filtered.
# 3. FILTERING STEPS:
#    a. Contaminants (CON_) and Reverse (REV_) hits were removed.
#    b. Only the top 500 proteins with highest intensity were kept. 
#       [O: "Only proteins present in all samples were kept"]
# 4. Columns were subsetted to include only: 'Protein IDs', 'Gene names', 
#    and 'LFQ intensity' columns.
# 5. The resulting file was saved and compressed using gzip.