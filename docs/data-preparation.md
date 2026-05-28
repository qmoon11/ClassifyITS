## Data Preparation

ClassifyITS does not run NCBI BLAST searches or build reference databases. These steps are already well implemented and should be run on HPC (high performance computing cluster). You must generate BLAST results externally and provide them as input. 
If preferred, equivalent R code using `system2()` is also provided after the shell workflow.
Below is a typical workflow to prepare your inputs (run on HPC):

1. Download BLAST+ (Linux example):
```bash
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.17.0+-x64-linux.tar.gz
    export PATH=/your/path/to/ncbi-blast-2.17.0+/bin:$PATH
```

2. Download the latest UNITE eukaryotic ITS database from <https://unite.ut.ee/>

    ClassifyITS does not download the UNITE database automatically. Users should download the appropriate UNITE FASTA file manually from the UNITE website before running BLAST.

    It is recommended to use the latest UNITE eukaryotic ITS database for best results. On the UNITE website, choose the "General FASTA release" that matches your analysis needs. For fungal ITS classification, the eukaryotic ITS release is typically appropriate.

    For example, you may download a file such as:

    ```text
    UNITE/euk/sh_general_release_dynamic_s_all_19.02.2025.fasta
    ```

    The exact filename may differ depending on the UNITE release date and selected database type. After downloading, place the FASTA file in your project directory or update the paths in the commands below to point to the downloaded file.

    For example, if your downloaded file is located at:

    ```text
    UNITE/euk/sh_general_release_dynamic_s_all_19.02.2025.fasta
    ```

    then use that path in the `makeblastdb` command below.
    
3. Build the BLAST database:
```bash
    makeblastdb \
      -in UNITE/euk/sh_general_release_dynamic_s_all_19.02.2025.fasta \
      -dbtype nucl \
      -out BLAST/uniteITSEuk_db
```
4. Run MegaBLAST search (example parameters):
```bash
    blastn \
      -task megablast \
      -query rep-seqs/dna-sequences.fasta \
      -db BLAST/uniteITSEuk_db \
      -out BLAST/megablast_ITS2.tsv \
      -word_size 28 \
      -reward 1 \
      -penalty -2 \
      -gapopen 0 \
      -gapextend 2 \
      -max_target_seqs 10 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
      -num_threads 8
```
This will generate a tab-delimited BLAST output file (for example, BLAST/megablast_ITS2.tsv) that can be used as input for ClassifyITS. It will contain the top 10 hits for each query sequence, along with relevant information such as percent identity, alignment length, e-value, and taxonomic information from the UNITE database.

**Required BLAST output format:**  
Make sure you use:  
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

Do not change the column order—this is required for ClassifyITS to work!

## Option 2: Prepare inputs from R

If you prefer to run the same preparation steps from R, you can use system2() to call the required Linux/BLAST commands. This does not mean ClassifyITS runs BLAST internally; rather, this is an example of how users can prepare the required BLAST output file themselves from an R session or R script.
For large datasets, this should still be run on an HPC compute node, for example through an Rscript submitted with a job scheduler such as SLURM.


# ============================================================
# Prepare BLAST output for ClassifyITS from R
# ============================================================
```r
# Number of CPU threads for BLAST
num_threads <- 8

# Input files
unite_fasta <- "UNITE/euk/sh_general_release_dynamic_s_all_19.02.2025.fasta"
query_fasta <- "rep-seqs/dna-sequences.fasta"

# Output directory and files
blast_dir <- "BLAST"
dir.create(blast_dir, recursive = TRUE, showWarnings = FALSE)

db_prefix <- file.path(blast_dir, "uniteITSEuk_db")
blast_output <- file.path(blast_dir, "megablast_ITS2.tsv")

# If BLAST+ is already available on your PATH, these should be found automatically.
makeblastdb <- Sys.which("makeblastdb")
blastn <- Sys.which("blastn")

# Alternatively, if you downloaded BLAST+ manually, specify the BLAST bin directory:
# blast_bin <- "/your/path/to/ncbi-blast-2.17.0+/bin"
# makeblastdb <- file.path(blast_bin, "makeblastdb")
# blastn <- file.path(blast_bin, "blastn")

# Check that BLAST tools are available
if (!nzchar(makeblastdb) || !file.exists(makeblastdb)) {
  stop("makeblastdb was not found. Install/load NCBI BLAST+ first.")
}

if (!nzchar(blastn) || !file.exists(blastn)) {
  stop("blastn was not found. Install/load NCBI BLAST+ first.")
}

# Check that input files exist
if (!file.exists(unite_fasta)) {
  stop("UNITE FASTA file was not found: ", unite_fasta)
}

if (!file.exists(query_fasta)) {
  stop("Query FASTA file was not found: ", query_fasta)
}

# Build the BLAST database
system2(
  makeblastdb,
  args = c(
    "-in", unite_fasta,
    "-dbtype", "nucl",
    "-out", db_prefix
  )
)

# Required BLAST output format for ClassifyITS
blast_outfmt <- paste(
  "6",
  "qseqid",
  "sseqid",
  "pident",
  "length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore",
  "stitle"
)

# Run MegaBLAST
system2(
  blastn,
  args = c(
    "-task", "megablast",
    "-query", query_fasta,
    "-db", db_prefix,
    "-out", blast_output,
    "-word_size", "28",
    "-reward", "1",
    "-penalty", "-2",
    "-gapopen", "0",
    "-gapextend", "2",
    "-max_target_seqs", "10",
    "-outfmt", blast_outfmt,
    "-num_threads", as.character(num_threads)
  )
)

message("BLAST output written to: ", blast_output)
```


If BLAST+ is not already available on your system, it can be downloaded manually before running the R workflow. For example, on Linux:

```r
dir.create("tools", recursive = TRUE, showWarnings = FALSE)

blast_url <- "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz"
blast_tar <- file.path("tools", basename(blast_url))

download.file(blast_url, blast_tar, mode = "wb")
untar(blast_tar, exdir = "tools")

blast_bin <- file.path("tools", "ncbi-blast-2.17.0+", "bin")

Sys.setenv(
  PATH = paste(blast_bin, Sys.getenv("PATH"), sep = .Platform$path.sep)
)

Sys.which("blastn")
Sys.which("makeblastdb")
```

Required BLAST output format
Make sure you use exactly:
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
Do not change the column order—this is required for ClassifyITS to work.