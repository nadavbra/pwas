from setuptools import setup
        
def readme():
    with open('README.rst', 'r') as f:
        return f.read()

setup(
    name = 'pwas',
    version = '1.0.4',
    description = 'Proteome-Wide Association Study (PWAS) is a protein-centric, gene-based method for conducting genetic association studies.',
    long_description = readme(),
    long_description_content_type = 'text/x-rst',
    url = 'https://github.com/nadavbra/pwas',
    author = 'Nadav Brandes',
    author_email  ='nadav.brandes@mail.huji.ac.il',
    packages = ['pwas', 'pwas.shared_utils'],
    license = 'MIT',
    scripts = [
        'bin/list_all_variants',
        'bin/determine_ref_alleles',
        'bin/organize_variant_effects_per_gene',
        'bin/calc_gene_effect_scores',
        'bin/pwas_test_genes',
        'bin/combine_pwas_results',
    ],
    install_requires = [
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'biopython',
        'statsmodels',
    ],
)
