## Optional: Custom Taxonomy Cutoffs

Optionally, you can provide a custom taxonomy cutoffs file via the `cutoffs_file` parameter in the pipeline.

Format:
The file should be a CSV or TSV table with columns for each rank, name and its cutoffs (either e-value or percent identity), for example:

taxonomic_rank,evalue_cutoff,pident_cutoff
kingdom,1e-5,80
phylum,1e-10,85

How to use:
- Specify the path to your custom cutoffs file using `cutoffs_file` in your pipeline call:

    ITS_assignment(
        blast_file = "input/blast_results.tsv",
        rep_fasta = "input/representative_seqs.fasta",
        cutoffs_file = "input/taxonomy_cutoffs.csv"
    )

- If omitted, the package default cutoffs will be used.

Example file:
A sample file is provided at `inst/extdata/taxonomy_cutoffs.csv` in package source code.

Reference example (in CSV format):

    taxonomic_rank,evalue_cutoff,pident_cutoff
    kingdom,1e-5,80
    phylum,1e-10,85
    class,1e-20,90
    order,1e-30,92
    family,1e-50,94
    genus,1e-70,96
    species,1e-80,98

Tip:
Check out `inst/extdata/taxonomy_cutoffs.csv` in the package folder for a full example and recommended formatting.