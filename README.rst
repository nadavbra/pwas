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

2. Genotype files [Currently the PLINK/BED and BGEN formats are supported. Work to support VCF is currently underway. It should be relatively easy to extend the code to support other formats as well.]


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
   
   
Step 1: Obtain the input genotype & phenotype files
---------------------------------------------------
