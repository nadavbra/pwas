import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import statsmodels.api as sm

from .shared_utils.util import log, normalize

# We stick with the default number of iterations (up to 35 iterations), but don't want to print any output (disp = 0)
FIT_PARAMS = dict(maxiter = 35, disp = 0)

class GeneTest(object):

    def __init__(self, phenotype_values, covariates, dominant_scores, recessive_scores, score_test_class):
        '''
        @param phenotype_values (pd.Series): A series of binary or continuous phenotypes (the response variable).
        @param covariates (pd.DataFrame): Covariates.
        @param dominant_scores & recessive_scores (np.array or pd.Series): dominant and recessive aggregated scores of the tested gene.
        @param score_test_class (BinaryTraitGeneScoreTest or ContinuousTraitGeneScoreTest): Depnding on whether the phenotype is binary or continuous.
        Note: phenotype_values, covariates, dominant_scores and recessive_scores must have compatible indices.
        '''
        self.phenotype_values = phenotype_values
        self.variables = covariates
        self.dominant_scores = dominant_scores
        self.recessive_scores = recessive_scores
        self.score_test_class = score_test_class
        
    def run(self, individually = True, combined = True):
        
        results = []
        
        if individually:
            results.append(self.run_individually())
        
        if combined:
            results.append(pd.Series([self.run_combined()], index = ['combined_pval']))
        
        return pd.concat(results)
    
    def run_individually(self):
        
        results = []
            
        for score_name, gene_scores in [('dominant', self.dominant_scores), ('recessive', self.recessive_scores)]:
            score_test = self.score_test_class(self.phenotype_values, gene_scores, self.variables)
            score_test.run()
            score_results = score_test.get_results()
            score_results.index = ['%s_%s' % (score_name, header) for header in score_results.index]
            results.append(score_results)
        
        return pd.concat(results)
        
    def run_combined(self):
    
        if self.dominant_scores.std() > 0 and self.recessive_scores.std() > 0:
            
            with pd.option_context('mode.chained_assignment', None):
                self.variables['gene_dominant_score'] = normalize(self.dominant_scores)
                self.variables['gene_recessive_scores'] = normalize(self.recessive_scores)
            
            combined_model = self.score_test_class.MODEL_CLASS(self.phenotype_values, self.variables)
            
            try:
                
                combined_model_results = combined_model.fit(**FIT_PARAMS)
                
                if not hasattr(combined_model_results, 'mle_retvals') or combined_model_results.mle_retvals['converged']:
                    contrast_mask = pd.DataFrame(0, index = [0, 1], columns = self.variables.columns)
                    contrast_mask.loc[0, 'gene_dominant_score'] = 1
                    contrast_mask.loc[1, 'gene_recessive_scores'] = 1
                    return combined_model_results.f_test(contrast_mask).pvalue
                else:
                    log('Logistic regression model failed to converge.')
            except np.linalg.LinAlgError:
                log('Failed fitting the model due to (probably) singular matrix.')
            finally:
                del self.variables['gene_dominant_score'], self.variables['gene_recessive_scores']
    
        return np.nan
    
class GeneScoreTest(object):
    
    def __init__(self, phenotype_values, gene_scores, variables, model_class):
        self.phenotype_values = phenotype_values
        self.gene_scores = gene_scores
        self.variables = variables
        self.model_class = model_class
        self.score_statistics = pd.Series([self.gene_scores.astype(np.float64).mean(), self.gene_scores.astype(np.float64).std()], \
                index = ['avg_score', 'std_score'])            
        self.model_results = None
        self.pval = np.nan
        
    def run(self):
        if self.score_statistics['std_score'] > 0:
        
            with pd.option_context('mode.chained_assignment', None):
                # It was found that normalizeing the scores helped the stability of logistic regression (preventing convergence failure)
                self.variables['gene_score'] = normalize(self.gene_scores)
        
            try:
                model = self.model_class(self.phenotype_values, self.variables)
                self.model_results = model.fit(**FIT_PARAMS)
                self.pval = self.model_results.pvalues['gene_score']
            except np.linalg.LinAlgError:
                log('Failed fitting the model due to (probably) singular matrix.')
            finally:
                del self.variables['gene_score']
            
class BinaryTraitGeneScoreTest(GeneScoreTest):

    MODEL_CLASS = sm.Logit
    
    def __init__(self, *args):
        super(BinaryTraitGeneScoreTest, self).__init__(*args, model_class = BinaryTraitGeneScoreTest.MODEL_CLASS)        
            
    def run(self):
        
        super(BinaryTraitGeneScoreTest, self).run()
        
        if self.model_results is not None and not self.model_results.mle_retvals['converged']:
            log('Logistic regression model failed to converge.')
            self.pval = np.nan
        
    def get_results(self):
        
        control_scores = self.gene_scores[self.phenotype_values == 0].astype(np.float64)
        avg_score_controls = control_scores.mean()
        std_score_controls = control_scores.std()
        
        case_scores = self.gene_scores[self.phenotype_values == 1].astype(np.float64)
        avg_score_cases = case_scores.mean()
        std_score_cases = case_scores.std()
        
        std_pooled = np.sqrt((std_score_controls ** 2 + std_score_cases ** 2) / 2)
        cohen_d = np.nan if std_pooled == 0 else (avg_score_cases - avg_score_controls) / std_pooled
        
        return pd.concat([self.score_statistics, pd.Series([avg_score_controls, std_score_controls, avg_score_cases, \
                std_score_cases, cohen_d, self.pval], index = ['avg_score_controls', 'std_score_controls', \
                'avg_score_cases', 'std_score_cases', 'cohen_d', 'pval'])])
                         
class ContinuousTraitGeneScoreTest(GeneScoreTest):

    MODEL_CLASS = sm.OLS
                         
    def __init__(self, *args):
        super(ContinuousTraitGeneScoreTest, self).__init__(*args, model_class = ContinuousTraitGeneScoreTest.MODEL_CLASS)
        
    def get_results(self):
        r, _ = pearsonr(self.phenotype_values, self.gene_scores)               
        return pd.concat([self.score_statistics, pd.Series([r, self.pval], index = ['r', 'pval'])])
