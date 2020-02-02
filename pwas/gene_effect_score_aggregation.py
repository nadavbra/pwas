'''
Common parameters:
- variant_genotypes <numpy float array of shape (n_variants, n_samples, n_alleles [3])>: For every variant and in the gene
(indexed by the first dimension) and every sample (indexed by the second dimension), what are the three probabilities of
it being either: i) homozygous allele-1, ii) heterozygous, iii) homozygous allele-2. The three probabilities must be
non-negative summing up to 1.
- effect_scores: <numpy float array of shape (n_variants,)>: The effect score (in the range 0 to 1) of every relevant
variant in the gene.
- allele1_refs: <numpy bool array of shape (n_variants,)>: Whether allele-1 is the reference variant (True) or allele-2
(False) for every relevant variant in the gene.
'''

import numpy as np

def calc_gene_dominant_effect_scores(variant_genotypes, effect_scores, allele1_refs, u, p, dtype = np.float32):
    x = calc_xy(variant_genotypes, effect_scores, allele1_refs, u, also_y = False, dtype = dtype)
    return calc_D(x, p)
    
def calc_gene_recessive_effect_scores(variant_genotypes, effect_scores, allele1_refs, u, p, q, dtype = np.float32):
    x, y = calc_xy(variant_genotypes, effect_scores, allele1_refs, u, dtype = dtype)
    return calc_R(x, y, p, q)

def calc_xy(variant_genotypes, effect_scores, allele1_refs, u, also_y = True, dtype = np.float32):

    '''
    The shape of x, y: (n_variants, n_samples)
    '''

    variant_genotypes = variant_genotypes.astype(dtype)
    effect_scores = effect_scores.astype(dtype)
    
    p = np.where(allele1_refs.reshape(-1, 1, 1), variant_genotypes, variant_genotypes[:, :, ::-1])
    s = effect_scores.reshape(-1, 1)
    
    p0 = p[:, :, 0]
    p1 = p[:, :, 1]
    p2 = p[:, :, 2]
    
    x = p0 + p1 * s + p2 * (u * s + (1 - u) * np.square(s))
    
    if also_y:
        y = p1 * (1 - s) + p2 * ((1 - u) * 2 * s * (1 - s))
        return x, y
    else:
        return x
    
def calc_D(x, p):
    minus_log_x = -np.log(x)
    lp_minus_log_x = np.linalg.norm(minus_log_x, p, axis = 0)
    return np.exp(-lp_minus_log_x)
    
def calc_R(x, y, p, q):
    
    x_zero_mask = (x == 0)
    num_zero_x = x_zero_mask.sum(axis = 0)
    no_zero_x_mask = (num_zero_x == 0)
    one_zero_x_mask = (num_zero_x == 1)
    
    R = np.zeros_like(num_zero_x, dtype = x.dtype)
    
    if one_zero_x_mask.any():
        R[one_zero_x_mask] = calc_R_with_one_zero_x(x[:, one_zero_x_mask], y[:, one_zero_x_mask], \
                x_zero_mask[:, one_zero_x_mask], p)
        
    if no_zero_x_mask.any():
        R[no_zero_x_mask] = calc_R_with_no_zero_x(x[:, no_zero_x_mask], y[:, no_zero_x_mask], p, q)
        
    return R

def calc_R_with_one_zero_x(x, y, x_zero_mask, p):
    if p == 1:
        return y[np.where(x_zero_mask)]
    else:
        return np.where(x_zero_mask, y, x).prod(axis = 0)
        
def calc_R_with_no_zero_x(x, y, p, q):
    zeta = np.linalg.norm(y / x, q, axis = 0)
    D = calc_D(x, p)
    return (zeta + 1) * D