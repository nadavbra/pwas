What is PWAS?
=============

Proteome-Wide Association Study (PWAS) is a protein-centric, gene-based method for conducting genetic association studies. PWAS detects protein-coding genes whose functional variabilities are correlated with given phenotypes across a cohort. It employs a machine-learning model to assess the functional damage caused to each protein within each sample (given the sample's genotype). These assessments are summarized as effect score matrices, where each combination of sample (row) and gene (column) is assigned a number between 0 (complete loss of function) to 1 (no effect). PWAS creates two such matrices, for either dominant or recessive inheritance. Following the creation of those matrices, PWAS can then test various phenotypes, looking for associations between the matrix columns (describing the functional variabilities of specific proteins) to phenotype values. In the case of a binary phenotype, a significant association would mean that the protein coded by the gene appears more damaged in cases than in controls (or vice versa).

For more details on PWAS you can refer to our paper `PWAS: proteome-wide association study—linking genes and phenotypes by functional variation in proteins <https://doi.org/10.1186/s13059-020-02089-x>`_, published in Genome Biology (2020).

Or, if you are more a video person, you can watch this talk on YouTube (originally given at RECOMB 2020):

.. image:: https://img.youtube.com/vi/TcgE_xb8ecw/0.jpg
   :target: https://www.youtube.com/watch?v=TcgE_xb8ecw


Installation
============

Dependencies
------------

PWAS requires Python 3.

Depending on the format of your genetic data, you will have to manually install a relevant Python parser before using PWAS. If you use the PLINK/BED format, you will have to install the `pandas-plink <https://pypi.org/project/pandas-plink/>`_ Python package. If you use the BGEN format, install `bgen_parser <https://github.com/nadavbra/bgen_parser>`_.

Part of PWAS's pipeline also requires other tools (see details below). Specifically, step 2.3 requires a variant assessment tool such as  `FIRM <https://github.com/nadavbra/firm>`_.

Upon installation, PWAS will automatically add the following Python packages:

* numpy
* scipy
* pandas
* matplotlib
* biopython
* statsmodels


Install PWAS
------------

Simply run:

.. code-block:: sh

   pip install pwas
   
**Important**: Make sure that the ``pip`` command refers to Python 3. If you are uncertain, consider using ``pip3`` instead.

Alternatively, to install PWAS directly from this GitHub repository, clone the project into a local directory and run from it:

.. code-block:: sh

   git submodule update --init
   python3 setup.py install


Obtaining the reference genome files
------------------------------------

If you need to determine the reference allele of each variant (step 2.2 described below), you will need to download the relevant version of the reference genome. Note that if you use FIRM you will have to download these files anyway. 

The reference genome sequences of all human chromosomes (chrXXX.fa.gz files) can be downloaded from UCSC's FTP site at: 

* ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ (for version hg19)
* ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ (for version hg38/GRCh38)

The chrXXX.fa.gz files need to be uncompressed to obtain chrXXX.fa files.

IMPORTANT: In version hg19 there's an inconsistency in the reference genome of the M chromosome between UCSC and RegSeq/GENCODE,
so the file chrM.fa has to be taken from RefSeq (NC_012920.1) instead of UCSC, from: https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&sort=&id=251831106&from=begin&to=end&maxplex=1. In GRCh38 all downloaded files should remain as they are.


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

   3.1. Collect the variant effect scores per gene
   
   3.2. Combine the variant effect scores with per-sample genotypes to obtain gene effect scores across the cohort's samples

4. **Find gene-phenotype associations**, which consists of:

   4.1. Run association tests (between a selected phenotype to the calculated gene effect scores)
   
   4.2. Collect the results and perform multiple-hypothesis testing correction
   
To ensure maximal flexibility and allow the integration of PWAS with other tools in a modular way, each of these steps consists of a separate command-line with well-defined inputs and outputs. This means that each of these steps can be skipped at your choice, given that you can provide the inputs necessary for the following steps by some alternative way.
   
   
Step 1: Obtain the input genotype & phenotype files
---------------------------------------------------

As stated, PWAS requires a CSV file with the phenotypic fields of your cohort. This CSV file requires a single column designated for unique sample identifiers (which should correspond to the identifiers in your genotype files). The CSV file should also contain one or more columns for the phenotypes you wish to test, and (preferably) covariates you wish to account for when testing the phenotypes (e.g. sex, age, genetic principal components, genetic batch, etc.). All phenotype and covariate fields must be numeric (i.e. 0s and 1s in the case of binary fields, or any number in the case of continuous fields).

If you work with the `UK Biobank <https://www.ukbiobank.ac.uk/>`_, you can use the `ukbb_parser package <https://github.com/nadavbra/ukbb_parser>`_ to easily create a CSV dataset with selected phenotype fields (and automatically extracted covariates for genetic association tests) through its `command-line interface <https://github.com/nadavbra/ukbb_parser#command-line-api>`_.

For example, the following command will create a suitable dataset with 49 prominent phenotypes (both binary/categorical and continuous) and 173 covariates extracted from the UK Biobank (assuming that you have access to the relevant UKBB fields).

.. code-block:: sh

    wget https://raw.githubusercontent.com/nadavbra/ukbb_parser/master/examples/phenotype_specs.py
    create_ukbb_phenotype_dataset --phenotype-specs-file=./phenotype_specs.py --output-dataset-file=./ukbb_dataset.csv --output-covariates-columns-file=./ukbb_covariate_columns.json

On top of the CSV of phenotypes, you will also need a CSV file specifying all the relevant genotyping files. This meta file is expected to list all the relevant genotype sources (one per row), having the following headers:

* **name**: A unique identifier of the genotype source (e.g. the name of the chromosome or genomic segment)
* **format**: The format of the genotype source (currently supporting only *plink* and *bgen*).

Genotype sources of *plink* format are expected to have three additional columns: **bed_file_path**, **bim_file_path** and **fam_file_path** (for the BED, BIM and FAM files, respectively). Likewise, genotype sources of *bgen* format are expected to have the following three columns: **bgen_file_path**, **bgi_file_path** and **sample_file_path** (for the .bgen, .bgen.bgi and .sample files, respectively).

Generating the meta CSV file of the genotype sources for the UK Biobank dataset can be easily achieved with the same ukbb_parser package. For example, the following command would generate the file for the imputated genotypes in BGEN format:

.. code-block:: sh

    create_ukbb_genotype_spec_file --genotyping-type=imputation --output-file=./ukbb_imputation_genotyping_spec.csv
    
**Very important note**: There's actually a good reason to choosing the UK Biobank's imputed genotypes over their raw markers. Unlike vanilla GWAS and other gene-based method (e.g. SKAT), for which it's sufficient to have some sampling of the variants in each Linkage Disequilibrium block, PWAS actually requires full knowledge of all the variants present in each sample. The underlying reason is that PWAS actually tries to figure out what happens to the genes (from functional perspective), and missing variants (with functional relevance) are likely to diminish its statistical power to uncover true associations. For this reason, PWAS is expected to work best with complete, unbiased genotyping (e.g. provided by whole-exome sequencing). If your genetic data was collected by SNP-array genotypes, then you will at least have to try to complete the missing variants through imputation.  
    
    
Step 2: Determine per-variant effect scores
-------------------------------------------


Step 2.1: List all the unique variants in the input genotyping files
--------------------------------------------------------------------

To combine all the variant descriptions across the input genotype sources into a unified list, simply use the ``list_all_variants`` command provided by PWAS.

For example, to list all the unique imputed variants in the UK Biobank, run:

.. code-block:: sh

    list_all_variants --genotyping-spec-file=./ukbb_imputation_genotyping_spec.csv --output-file=./ukbb_imputed_variants.csv --verbose


Step 2.2 (optional): Determine the reference allele of each variant
-------------------------------------------------------------------

In most genetic datasets it is the convention that the first allele listed in each variant is the reference allele and the second is the alternative alleles. However, in some datasets (including the UK Biobank) this convention is sometimes broken. In order to function properly, PWAS needs to know which of the two alleles listed in each variant is the reference allele. If you are not sure whether this convention holds in your dataset, it is recommended that you determine the reference alleles, just to be on the safe side. The ``determine_ref_alleles`` command (provided by PWAS) will compare each variant against the reference genome to validate which of the two variants is the reference allele.

For example, to determine the reference alleles of the imputed UKBB variants, run:

.. code-block:: sh

    determine_ref_alleles --variants-file=./ukbb_imputed_variants.csv --ref-genome-dir=/path/to/hg19/ --chrom-col=chromosome --pos-col=position --allele1-col=allele1 --allele2-col=allele2 --override --verbose
    
where the --ref-genome-dir option should point to a directory with the sequences of the relevant version of the human reference genome (hg19 in the case of the UKBB). This directory is expected to have one (uncompressed) FASTA file per chromosome (e.g. chr1.fa, chr2.fa, ..., chr22.fa, chrX.fa, chrY.fa, chrM.fa). See the `Obtaining the reference genome files <#obtaining-the-reference-genome-files>`_ section above.


Step 2.3: Calculate the effect score of each variant
----------------------------------------------------

A crucial step in determining the functional status of genes is to first determine the predicted functional effects of individual variants. PWAS requires that each variant will be assigned an effect score between 0 (indicating complete loss of function of the gene) to 1 (indicating no effect). PWAS has been designed and tested to work with `FIRM <https://github.com/nadavbra/firm>`_, a machine-learning framework for predicting the functional impact of variants affecting protein sequences at the molecular-level. However, PWAS is completely generic and could, in principle, work with any variant assessment tool (e.g. `CADD <https://cadd.gs.washington.edu/>`_). In fact, since all of PWAS's calculations are derived from the per-variant effect scores, and it's actually agnostic to their interpretation, you can even assign scores to non-coding genes or use scores that capture other biological properties of mutations (even though PWAS was originally designed for discovering proteomic associations).

Whatever tool you end up using, you will need to produce a `JSON-lines <http://jsonlines.org/>`_ file. Each row in the file is expected to describe the effects of the variants in the corresponding row in the variants CSV file (in particular, the two files are expected to have the same number of lines, except the headers line that is only expected in the CSV file, but not in the JSON-lines file). Each row in the file is expected to be a JSON-formatted dictionary, mapping each gene index (a running integer index arbitrarily assigned to each gene) into the variant's list of effects on the gene, each is a pair of i) effect description (string) and ii) effect score (float, between 0 to 1).

For example, to calculate the effect scores of UKBB's imputed variants with FIRM (following its installation), run:

.. code-block:: sh

    firm_determine_extended_gene_effects_and_scores --variants-csv-file=./ukbb_imputed_variants.csv --output-effects-file=./ukbb_imputation_effects.jsonl --genes-dir=./ --ref-genome=GRCh37 --chrom-col=chromosome --pos-col=position --allele1-col=allele1 --allele2-col=allele2 --is-allele1-ref-col=is_allele1_ref
    
    
Step 3: Aggregate the per-variant into per-gene effect scores
-------------------------------------------------------------


Step 3.1: Collect the variant effect scores per gene
----------------------------------------------------

Having completed step 2, you should now have: i) a CSV file listing all the variants genotyped in your cohort, and ii) a JSON-lines file specifying all the effects of these variants on genes, where each variant-gene effect is assigned a functional score. In order to aggregate the per-variant effect scores into per-gene scores, PWAS first needs the variant effects to be organized per gene. It requires a separate CSV file per gene listing all the variants affecting that gene. These CSV files should have, on top of all the columns in the original CSV file (that lists all the variants), an additional *effect_score* column with the effect score of each of the variants (with respect to the file's gene).

To generate the per-gene files, simply use the ``organize_variant_effects_per_gene`` command provided by PWAS.

For example, the following will generate the required per-gene CSV files for the imputed variants in the UKBB:

.. code-block:: sh

    mkdir ./ukbb_imputation_variants_per_gene
    organize_variant_effects_per_gene --variants-file=./ukbb_imputed_variants.csv --effects-file=./ukbb_imputation_effects.jsonl --gene-variants-dir=./ukbb_imputation_variants_per_gene/
    
    
Step 3.2: Calculate the gene effect scores
------------------------------------------

Now here comes PWAS's magic sauce. We are going to aggregate the per-variant effect scores into per-gene (dominant and recessive) effect scores, while taking into account each sample's unique genotype. The relevant PWAS command is ``calc_gene_effect_scores``.

For example, the following command will calculate the gene effect scores for all of the UK Biobank's samples, based on their imputed genotypes:

.. code-block:: sh

   mkdir ./ukbb_imputation_gene_effect_scores/
   calc_gene_effect_scores --genotyping-spec-file=./ukbb_imputation_genotyping_spec.csv --gene-variants-dir=./ukbb_imputation_variants_per_gene/ --gene-effect-scores-dir=./ukbb_imputation_gene_effect_scores/ --is-allele1-ref-col=is_allele1_ref

Since this process is computationally intensive (with respect to storage and CPU), it might be a good idea to distribute it across multiple tasks (and potentially sending them to run on a cluster). Luckily for you, this command is already equipped with built-in distribution functionality. For a full explanation on all the different options to distribute the command, please refer to its help message. 

In our example, we can distribute the process into 1,000 tasks and send them to run on a cluster managed by SLURM, by running:

.. code-block:: sh

   sbatch --array=0-999 --mem=32g -c1 --time=1-0 --wrap="calc_gene_effect_scores --genotyping-spec-file=./ukbb_imputation_genotyping_spec.csv --gene-variants-dir=./ukbb_imputation_variants_per_gene/ --gene-effect-scores-dir=./ukbb_imputation_gene_effect_scores/ --is-allele1-ref-col=is_allele1_ref --task-index-env-variable=SLURM_ARRAY_TASK_ID --total-tasks-env-variable=SLURM_ARRAY_TASK_COUNT"
   
Once the jobs have successfully finished, you should have a CSV file per gene, with the effect scores of each sample.

It might be a good idea to validate that you have the correct number of CSV files (i.e. the same as the number of CSV files listing the per-gene variants):

.. code-block:: sh

   ls -l ./ukbb_imputation_variants_per_gene/ | wc -l
   ls -l ./ukbb_imputation_gene_effect_scores/ | wc -l
   
The algorithm that aggregates the variant effect scores into gene effect scores is actually dependent on 5 parameters that the ``calc_gene_effect_scores`` command allows you to specify, although the default values are likely a sensible choice. For the full mathematical details of the aggregation algorithm, and the meaning of those parameters, please refer to our paper.


Step 4: Find gene-phenotype associations
----------------------------------------


Step 4.1: Run association tests
--------------------------------

Having gone through step 1, you should have a CSV file with phenotypes and covariates, and having completed step 3 you should also have per-gene CSV files with the gene effect scores. The last step of PWAS is to simply look for statistical correlations between the phenotypes to the gene scores, in order to uncover gene-phenotype associations (with respect to the functional variability captured by the pre-calculated gene effect scores, which, in the default case where FIRM has been used as the variant assessment tool, reflect the estimated functions of the proteins coded by those genes). In fact, this step consists of nothing more than routine statistical methods (linear and logistic regression), and you could, in principle, use any statistics software of your choice (e.g. PLINK, R, etc.). Still, PWAS comes with its own built-in implementation which also provides, on top p-values, some additional unique metrics. Unless you feel very confident that you know what you are doing, it is recommended that you just use the implementation of PWAS, as provided by the ``pwas_test_genes`` command.

To continue our ongoing UKBB example, let's say we want to find PWAS associations for type-II diabetes. Then simply run:

.. code-block:: sh

   mkdir ./ukbb_imputation_per_gene_type2_diabetes_pwas_results
   pwas_test_genes --dataset-file=./ukbb_dataset.csv --gene-effect-scores-dir=./ukbb_imputation_gene_effect_scores/ --per-gene-pwas-results-dir=./ukbb_imputation_per_gene_type2_diabetes_pwas_results/ --sample-id-col=eid --phenotype-col="Type 2 diabetes" --covariate-cols-json-file=./ukbb_covariate_columns.json
   
This process will go through each gene in ``./ukbb_imputation_gene_effect_scores/`` and run a logistic regression test of the "Type 2 diabetes" column in ``./ukbb_dataset.csv`` against the gene's effect scores (while also taking into account the covariates in the columns specified by ``./ukbb_covariate_columns.json``). It will save the resulted summary statistics of each gene as a separate CSV file in ``./ukbb_imputation_per_gene_type2_diabetes_pwas_results/``.

This process too can be computationally intensive (in terms of CPU time), especially for large datasets (with many samples and covariates) such as the UKBB. Fortunately, the ``pwas_test_genes`` command comes with a built-in functionality that allows one to distribute it across many computing resources. For full details on that, please refer to its help message. As an example, if you want to distribute the process across 1,000 tasks and send them to run on a cluster managed by SLURM, simply run:

.. code-block:: sh

   sbatch --array=0-999 --mem=32g -c1 --time=1-0 --wrap="pwas_test_genes --dataset-file=./ukbb_dataset.csv --gene-effect-scores-dir=./ukbb_imputation_gene_effect_scores/ --per-gene-pwas-results-dir=./ukbb_imputation_per_gene_type2_diabetes_pwas_results/ --sample-id-col=eid --phenotype-col='Type 2 diabetes' --covariate-cols-json-file=./ukbb_covariate_columns.json --task-index-env-variable=SLURM_ARRAY_TASK_ID --total-tasks-env-variable=SLURM_ARRAY_TASK_COUNT"
   
Here too, once everything is done and over with, it will be a good idea to validate that you've got the right number of files. These two command are expected to give you the same number:
   
.. code-block:: sh

   ls -l ./ukbb_imputation_gene_effect_scores/ | wc -l
   ls -l ./ukbb_imputation_per_gene_type2_diabetes_pwas_results/ | wc -l
   
   
Step 4.2: Collect the results and perform multiple-hypothesis testing correction
--------------------------------------------------------------------------------

To collect the summary statistics calculated in the previous step (which are currently spread across many CSV files) and perform multiple-hypothesis testing correction, simply use the ``combine_pwas_results`` command.

In our ongoing example, just run:

.. code-block:: sh

   combine_pwas_results --genes-file=./genes_hg19.csv --per-gene-pwas-results-dir=./ukbb_imputation_per_gene_type2_diabetes_pwas_results/ --results-file=./ukbb_imputation_type2_diabetes_pwas_results.csv
   
The file ``./genes_hg19.csv`` should have been generated by FIRM when you used it to estimate the variant effect scores. It is necessary to provide the details of all the genes which, up until this point, PWAS represented by nothing more than indices. PWAS is actually agnostic to the content of this file, and it simply concatenates it, as is, before the summary statistics of each gene.

When the process is finished, you will have the file ``./ukbb_imputation_type2_diabetes_pwas_results.csv`` with the complete summary statistics of all tested genes. And that's the end of it - you are now the proud owner of freshly generated PWAS results!


Recap (the complete pipeline)
-----------------------------

For quicker future reference, here's the complete pipeline for running PWAS for type-II diabetes over the imputed genotypes provided by the UK Biobank.

First, to generate the necessary phenotype & genotype files from the UKBB dataset, use the ``ukbb_parser`` package:

.. code-block:: sh

   wget https://raw.githubusercontent.com/nadavbra/ukbb_parser/master/examples/phenotype_specs.py
   create_ukbb_phenotype_dataset --phenotype-specs-file=./phenotype_specs.py --output-dataset-file=./ukbb_dataset.csv --output-covariates-columns-file=./ukbb_covariate_columns.json
   create_ukbb_genotype_spec_file --genotyping-type=imputation --output-file=./ukbb_imputation_genotyping_spec.csv
   
Second, you will have to list all the dataset's variants and determine the reference allele of each variant:

.. code-block:: sh

   list_all_variants --genotyping-spec-file=./ukbb_imputation_genotyping_spec.csv --output-file=./ukbb_imputed_variants.csv --verbose
   determine_ref_alleles --variants-file=./ukbb_imputed_variants.csv --ref-genome-dir=/path/to/hg19/ --chrom-col=chromosome --pos-col=position --allele1-col=allele1 --allele2-col=allele2 --override --verbose
   
And then calculate the variant effect scores (here using FIRM):

.. code-block:: sh

   firm_determine_extended_gene_effects_and_scores --variants-csv-file=./ukbb_imputed_variants.csv --output-effects-file=./ukbb_imputation_effects.jsonl --genes-dir=./ --ref-genome=GRCh37 --chrom-col=chromosome --pos-col=position --allele1-col=allele1 --allele2-col=allele2 --is-allele1-ref-col=is_allele1_ref
   
Next, you will need to organize the variant effect scores per gene and aggregate them into gene effect scores (distributing that process on a cluster to speed things up):

.. code-block:: sh

   mkdir ./ukbb_imputation_variants_per_gene
   organize_variant_effects_per_gene --variants-file=./ukbb_imputed_variants.csv --effects-file=./ukbb_imputation_effects.jsonl --gene-variants-dir=./ukbb_imputation_variants_per_gene/
   mkdir ./ukbb_imputation_gene_effect_scores/
   sbatch --array=0-999 --mem=32g -c1 --time=1-0 --wrap="calc_gene_effect_scores --genotyping-spec-file=./ukbb_imputation_genotyping_spec.csv --gene-variants-dir=./ukbb_imputation_variants_per_gene/ --gene-effect-scores-dir=./ukbb_imputation_gene_effect_scores/ --is-allele1-ref-col=is_allele1_ref --task-index-env-variable=SLURM_ARRAY_TASK_ID --total-tasks-env-variable=SLURM_ARRAY_TASK_COUNT"
   
And validate that you got the correct number of files:

.. code-block:: sh

   ls -l ./ukbb_imputation_variants_per_gene/ | wc -l
   ls -l ./ukbb_imputation_gene_effect_scores/ | wc -l
   
Lastly, run the actual association tests (again using a cluster):

.. code-block:: sh

   mkdir ./ukbb_imputation_per_gene_type2_diabetes_pwas_results
   sbatch --array=0-999 --mem=32g -c1 --time=1-0 --wrap="pwas_test_genes --dataset-file=./ukbb_dataset.csv --gene-effect-scores-dir=./ukbb_imputation_gene_effect_scores/ --per-gene-pwas-results-dir=./ukbb_imputation_per_gene_type2_diabetes_pwas_results/ --sample-id-col=eid --phenotype-col='Type 2 diabetes' --covariate-cols-json-file=./ukbb_covariate_columns.json --task-index-env-variable=SLURM_ARRAY_TASK_ID --total-tasks-env-variable=SLURM_ARRAY_TASK_COUNT"
   
Validate that you've got all the files:

.. code-block:: sh
   
   ls -l ./ukbb_imputation_gene_effect_scores/ | wc -l
   ls -l ./ukbb_imputation_per_gene_type2_diabetes_pwas_results/ | wc -l
    
And combine the results to get the final summary statistics file:

.. code-block:: sh
   
   combine_pwas_results --genes-file=./genes_hg19.csv --per-gene-pwas-results-dir=./ukbb_imputation_per_gene_type2_diabetes_pwas_results/ --results-file=./ukbb_imputation_type2_diabetes_pwas_results.csv
 

License
=======
PWAS is a free open source project available under the `MIT License <https://en.wikipedia.org/wiki/MIT_License>`_.
 
   
Cite us
=======

If you use PWAS as part of work contributing to a scientific publication, we ask that you cite our paper: Brandes, N., Linial, N. & Linial, M. PWAS: proteome-wide association study—linking genes and phenotypes by functional variation in proteins. Genome Biol 21, 173 (2020). https://doi.org/10.1186/s13059-020-02089-x
