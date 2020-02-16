What is PWAS?
=============

Proteome-Wide Association Study (PWAS) is a protein-centric, gene-based method for conducting genetic association studies. PWAS detects protein-coding genes whose functional variabilities are correlated with given phenotypes across a cohort. It employs a machine-learning model to assess the functional damage caused to each protein within each sample (given the sample's genotype). These assessments are summarized as effect score matrices, where each combination of sample (row) and gene (column) is assigned a number between 0 (complete loss of function) to 1 (no effect). PWAS creates two such matrices, for either dominant or recessive inharitance. Following the creation of those matrices, PWAS can then test various phenotypes, looking for associations between the matrix columns (describing the functional variabilities of specific proteins) to phenotype values. In the case of a binary phenotype, a significant association would mean that the protein coded by the gene appears more damaged in cases than in controls (or vice versa).

For more details read our manuscript: Nadav Brandes, Nathan Linial, Michal Linial, PWAS: Proteome-Wide Association Study, bioRxiv, https://doi.org/10.1101/812289


Usage
=====

Overview
--------


PWAS requires the following input files:

1. Phenotypes and (optionally) covariates in a CSV file

2. Genotype files [Currently only the PLINK/BED and BGEN formats are supported. An effort to also support VCF files is currently underway, and it should be relatively easy to extend the code to support other formats as well.]


Running PWAS consists of the following steps:

1. **Obtain the input genotype & phenotype files**

2. **Determine per-variant effect scores**, which consists of:

   2.1. List all the unique variants in the input genotyping files
  
   2.2. (Optional) Determine the reference allele of each variant
  
   2.3. Calculate the effect score of each variant (using the variant assessment tool of your choice)

3. **Aggregate the per-variant into per-gene effect scores**, which consists of:

   3.1. Collect the varaint effect scores per gene
   
   3.2. Combine the variant effect scores with per-sample genotypes to obtain gene effect scores across the cohort's samples

4. **Run the association tests**, which consists of:

   4.1. Run the statistical associations tests (between a selected phenotype to the calculated gene effect scores)
   
   4.2. Collect the results and perform multiple-hypothesis testing correction
   
To ensure maximal flexibility and allow the integration of PWAS with other tools in a modular way, each of these steps consists of a separate command-line with well-defined inputs and outputs. This means that each of these steps can be skipped at your choice, given that you can provide the inputs necessary for the following steps by some alternative way.
   
   
Step 1: Obtain the input genotype & phenotype files
---------------------------------------------------

As stated, PWAS requires a CSV file with the phenotypic fields of your cohort. This CSV file requires a single column designated for unique sample identifiers (which should correspond to the identifiers in your genotype files). The CSV file should also contain one or more columns for the phenotypes you wish to test, and (preferably) covariates you wish to account for when testing the phenotypes (e.g. sex, age, genetic principal components, genetic batch, etc.). All phenotype and covariate fields must be numeric (i.e. 0s and 1s in the case of binary fields, or any number in the case of continuous fields).

If you work with the `UK Biobank <https://www.ukbiobank.ac.uk/>`_, you can use the `ukbb_parser package <https://github.com/nadavbra/ukbb_parser>`_ to easily create a CSV dataset with selected phenotype fields (and automatically extracted covariates for genetic association tests) through its `command-line interface <https://github.com/nadavbra/ukbb_parser#command-line-api>`_.

For example, the following command will create a suitable dataset with 49 prominent phenotypes (both binary/categorical and continuous) and 173 covariates extracted from the UK Biobank (assuming that you have access to the relevant UKBB fields).

.. code-block:: cshell

    wget https://raw.githubusercontent.com/nadavbra/ukbb_parser/master/examples/phenotype_specs.py
    create_ukbb_phenotype_dataset --phenotype-specs-file=./phenotype_specs.py --output-dataset-file=./ukbb_dataset.csv --output-covariates-columns-file=./ukbb_covariate_columns.json

On top of the CSV of phenotypes, you will also need a CSV file specifying all the relevant genotyping files. This meta file is expected to list all the relevant genotype sources (one per row), having the following headers:

* **name**: A unique identifier of the genotype source (e.g. the name of the chromosome or genomic segment)
* **format**: The format of the genotype source (currently supporting only *plink* and *bgen*).

Genotype sources of *plink* format are expected to have three additional columns: **bed_file_path**, **bim_file_path** and **fam_file_path** (for the BED, BIM and FAM files, respectively). Likewise, genotype sources of *bgen* format are expected to have the following three columns: **bgen_file_path**, **bgi_file_path** and **sample_file_path** (for the .bgen, .bgen.bgi and .sample files, respectively).

Generating the meta CSV file of the genotype sources for the UK Biobank dataset can be easily achieved with the same ukbb_parser package. For example, the following command would generate the file for the imputated genotypes in BGEN format:

.. code-block:: cshell

    create_ukbb_genotype_spec_file --genotyping-type=imputation --output-file=./ukbb_imputation_genotyping_spec.csv
    
    
Step 2: Determine per-variant effect scores
-------------------------------------------


Step 2.1: List all the unique variants in the input genotyping files
--------------------------------------------------------------------

To combine all the varaint descriptions across the input genotype sources into a unified list, simply use the ``list_all_variants`` command provided by PWAS.

For example, to list all the unique imputed variants in the UK Biobank, run:

.. code-block:: cshell

    list_all_variants --genotyping-spec-file=./ukbb_imputation_genotyping_spec.csv --output-file=./ukbb_imputed_variants.csv --verbose


Step 2.2 (optional): Determine the reference allele of each variant
-------------------------------------------------------------------

In most genetic datasets it is the convention that the first allele listed in each variant is the reference allele and the second is the alternative alleles. However, in some datasets (including the UK Biobank) this convention is sometimes broken. In order to function properly, PWAS needs to know which of the two alleles listed in each variant is the reference allele. If you are not sure whether this convention holds in your dataset, it is recommended that you determine the reference alleles, just to be on the safe side. The ``determine_ref_alleles`` command (provided by PWAS) will compare each variant against the reference genome to validate which of the two variants is the reference allele.

For example, to determine the reference alleles of the imputed UKBB variants, run:

.. code-block:: cshell

    determine_ref_alleles --variants-file=./ukbb_imputed_variants.csv --ref-genome-dir=/path/to/hg19/ --chrom-col=chromosome --pos-col=position --allele1-col=allele1 --allele2-col=allele2 --override --verbose
    
where the --ref-genome-dir option should point to a directory with the sequences of the relevant version of the human reference genome (hg19 in the case of the UKBB). This directory is expected to have one (uncompressed) FASTA file per chromosome (e.g. chr1.fa, chr2.fa, ..., chr22.fa, chrX.fa, chrY.fa, chrM.fa). See the `Obtaining the reference genome files <#obtaining-the-reference-genome-files>`_ section below.


Step 2.3: Calculate the effect score of each variant
----------------------------------------------------

A crucial step in determining the functional status of genes is to first determine the predicted functional effects of individual variants. PWAS requires that each variant will be assigned an effect score between 0 (indicating complete loss of function of the gene) to 1 (indicating no effect). PWAS has been designed and tested to work with `FIRM <https://github.com/nadavbra/firm>`_, a machine-learning framework for predicting the functional impact of variants affecting protein sequences at the molecular-level. However, PWAS is completely generic and could, in principle, work with any variant assessment tool (e.g. `CADD <https://cadd.gs.washington.edu/>`_). In fact, since all of PWAS's calculations are derived from the per-variant effect scores, and it's actually agnostic to their interpretation, you can even assign scores to non-coding genes or use scores that capture other biological properties of mutations (even though PWAS was originally designed for discovering proteomic associations).

Whatever tool you end up using, you will need to produce a `JSON-lines <http://jsonlines.org/>`_ file. Each row in the file is expected to describe the effects of the variants in the corresponding row in the variants CSV file (in particular, the two files are expected to have the same number of lines, except the headers line that is only expected in the CSV file, but not in the JSON-lines file). Each row in the file is expected to be a JSON-formatted dictionary, mapping each gene index (a running integer index arbitrarily assigned to each gene) into the variant's list of effects on the gene, each is a pair of i) effect description (string) and ii) effect score (float, between 0 to 1).

For example, to calculate the effect scores of UKBB's imputed variants with FIRM (following its installation), run:

.. code-block:: cshell

    firm_determine_extended_gene_effects_and_scores --variants-csv-file=./ukbb_imputed_variants.csv --output-effects-file=./ukbb_imputation_effects.jsonl --genes-dir=./ --ref-genome=GRCh37 --chrom-col=chromosome --pos-col=position --allele1-col=allele1 --allele2-col=allele2 --is-allele1-ref-col=is_allele1_ref
    
    
Step 3: Aggregate the per-variant into per-gene effect scores
-------------------------------------------------------------


Step 3.1: Collect the varaint effect scores per gene
----------------------------------------------------

Having completed step 2, you should now have: i) a CSV file listing all the variants genotyped in your cohort, and ii) a JSON-lines file specifying all the effects of these variants on genes, where each variant-gene effect is assigned a functional score. In order to aggregate the per-variant effect scores into per-gene scores, PWAS first needs the variant effects to be organized per gene. It requires a seperate CSV file per gene listing all the variants affecting that gene. These CSV files should have, on top of all the columns in the original CSV file (that lists all the variants), an additional *effect_score* column with the effect score of each of the variants (with respect to the file's gene).

To generate the per-gene files, simply use the ``organize_variant_effects_per_gene`` command provided by PWAS.

For example, the following will generate the required per-gene CSV files for the imputed variants in the UKBB:

.. code-block:: cshell

    mkdir ./ukbb_imputation_variants_per_gene
    organize_variant_effects_per_gene --variants-file=./ukbb_imputed_variants.csv --effects-file=./ukbb_imputation_effects.jsonl --gene-variants-dir=./ukbb_imputation_variants_per_gene/
    
    
Step 3.2: Calculate the gene effect scores
------------------------------------------

Now here comes PWAS's magic sauce. We are going to aggregate the per-variant effect scores into per-gene (dominant and recessive) effect scores, while taking into account each sample's unique genotype. The relevant PWAS command is ``calc_gene_effect_scores``.

For example, the following command will calculate the gene effect scores for all of the UK Biobank's samples, based on their imputed genotypes:

.. code-block:: cshell

   mkdir ./ukbb_imputation_gene_effect_scores/
   calc_gene_effect_scores --genotyping-spec-file=./ukbb_imputation_genotyping_spec.csv --gene-variants-dir=./ukbb_imputation_variants_per_gene/ --gene-effect-scores-dir=./ukbb_imputation_gene_effect_scores/ --is-allele1-ref-col=is_allele1_ref

Since this process is computationally intensive (with respect to storage and CPU), it might be a good idea to distribute it across multiple tasks (and potentially sending them to run on a cluster). Luckily for you, this command is already equipped with built-in distribution functionality. For a full explanation on all the different options to distribute the command, please refer to its help message. 

In our example, we can distribute the process into 1,000 tasks and send them to run on a cluster managed by SLURM, by running:

.. code-block:: cshell

   sbatch --array=0-999 --mem=32g -c1 --time=1-0 --wrap="calc_gene_effect_scores --genotyping-spec-file=./ukbb_imputation_genotyping_spec.csv --gene-variants-dir=./ukbb_imputation_variants_per_gene/ --gene-effect-scores-dir=./ukbb_imputation_gene_effect_scores/ --is-allele1-ref-col=is_allele1_ref --task-index-env-variable=SLURM_ARRAY_TASK_ID --total-tasks-env-variable=SLURM_ARRAY_TASK_COUNT"
   
Once the jobs have successfully finished, you should have a CSV file per gene, with the effect scores of each sample.

It might be a good idea to validate that you have the correct number of CSV files (i.e. the same as the number of CSV files listing the per-gene variants):

.. code-block:: cshell

   ls -l ./ukbb_imputation_variants_per_gene/ | wc -l
   ls -l ./ukbb_imputation_gene_effect_scores/ | wc -l
   
The algorithm that aggregates the variant effect scores into gene effect scores is actually dependent on 5 parameters that the ``calc_gene_effect_scores`` command allows you to specifiy, although the default values are likely a sensible choice. For the full mathematical details of the aggregation algorithm, and the meaning of those parameters, please refer to our paper.


Installation
============

Obtaining the reference genome files
------------------------------------

The reference genome sequences of all human chromosomes (chrXXX.fa.gz files) can be downloaded from UCSC's FTP site at: 

* ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ (for version hg19)
* ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ (for version hg38/GRCh38)

The chrXXX.fa.gz files need to be uncompressed to obtain chrXXX.fa files.

IMPORTANT: In version hg19 there's an inconsistency in the reference genome of the M chromosome between UCSC and RegSeq/GENCODE,
so the file chrM.fa has to be taken from RefSeq (NC_012920.1) instead of UCSC, from: https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&sort=&id=251831106&from=begin&to=end&maxplex=1. In GRCh38 all downloaded files should remain as they are.




