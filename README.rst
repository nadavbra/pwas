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


