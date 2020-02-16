What is PWAS?
=============

Proteome-Wide Association Study (PWAS) is a protein-centric, gene-based method for conducting genetic association studies. PWAS detects protein-coding genes whose functional variabilities are correlated with given phenotypes across a cohort. It employs a machine-learning model to assess the functional damage caused to each protein within each sample (given the sample's genotype). These assessments are summarized as effect score matrices, where each combination of sample (row) and gene (column) is assigned a number between 0 (complete loss of function) to 1 (no effect). PWAS creates two such matrices, for either dominant or recessive inharitance. Following the creation of those matrices, PWAS can then test various phenotypes, looking for associations between the matrix columns (describing the functional variabilities of specific proteins) to phenotype values. In the case of a binary phenotype, a significant association would mean that the protein coded by the gene appears more damaged in cases than in controls (or vice versa).

For more details read our manuscript: Nadav Brandes, Nathan Linial, Michal Linial, PWAS: Proteome-Wide Association Study, bioRxiv, https://doi.org/10.1101/812289


Usage
=====
