import os
import inspect
from getpass import getuser

import numpy as np
import pandas as pd

from .shared_utils.util import safe_symlink

class GenotypingManager:
    
    def __init__(self, genotyping_spec_df, temp_dir):
        self.genotyping_spec_df = genotyping_spec_df
        self.temp_dir = temp_dir
        self._index_to_reader = {}
        
    def get_genotyping_reader(self, index):
        
        if index not in self._index_to_reader:
            self._index_to_reader[index] = _get_reader_for_genotyping_spec_record(self.genotyping_spec_df.loc[index], self.temp_dir)
        
        return self._index_to_reader[index]
        
    def get_genotyping_reader_by_name(self, name):
        index, = self.genotyping_spec_df[self.genotyping_spec_df['name'] == name].index
        return self.get_genotyping_reader(index)
        
    def iter_genotyping_readers(self):
        for index, genotyping_spec_record in self.genotyping_spec_df.iterrows():
            yield index, genotyping_spec_record, self.get_genotyping_reader(index)
            
    def __len__(self):
        return len(self.genotyping_spec_df)

class BgenGenotypingReader:

    def __init__(self, bgen_file_path, bgi_file_path, sample_file_path):
    
        try:
            from bgen_parser import BgenParser
        except ImportError:
            raise ImportError('Failed importing bgen_parser.BgenParser. Make sure bgen_parser is installed. See: https://github.com/nadavbra/bgen_parser.')

        self.bgen_parser = BgenParser(bgen_file_path, bgi_file_path, sample_file_path)
        
    def get_variants(self):
        return self.bgen_parser.variants
        
    def get_sample_ids(self):
        return self.bgen_parser.sample_ids
        
    def get_variant_probs(self, variant_index):
        return self.bgen_parser.read_variant_probs(variant_index)
        
class PlinkGenotypingReader:

    def __init__(self, bed_file_path, bim_file_path, fam_file_path, temp_dir):
        
        try:
            from pandas_plink import read_plink
        except ImportError:
            raise ImportError('Failed importing pandas_plink.read_plink. Make sure pandas-plink is installed. See: https://pypi.org/project/pandas-plink/.')
            
        plink_path_prefix = _create_plink_links(bed_file_path, bim_file_path, fam_file_path, temp_dir)
        self.bim, self.fam, self.G = read_plink(plink_path_prefix)
        
    def get_variants(self):
        return self.bim
        
    def get_sample_ids(self):
        return self.fam['iid']
        
    def get_variant_probs(self, variant_index):
        variant_genotypes = self.G[variant_index, :].compute()
        assert set(np.unique(variant_genotypes[~np.isnan(variant_genotypes)])) <= {0, 1, 2}
        return np.arange(3).reshape(1, 3) == variant_genotypes.reshape(-1, 1)

def _get_reader_for_genotyping_spec_record(genotyping_spec_record, temp_dir):
    reader_type = _GENOTYPING_FORMAT_TO_READER_TYPE[genotyping_spec_record['format']]
    all_reader_args = set(inspect.getfullargspec(reader_type.__init__).args[1:])
    optional_reader_arg_assignments_from_other_sources = dict(temp_dir = temp_dir)
    reader_args_from_spec = {arg for arg in all_reader_args if arg not in optional_reader_arg_assignments_from_other_sources}
    reader_arg_assignments_from_spec = genotyping_spec_record[list(reader_args_from_spec)].to_dict()
    rader_arg_assignments_from_other_sources = {arg: value for arg, value in optional_reader_arg_assignments_from_other_sources.items() if \
            arg in all_reader_args}
    return reader_type(**reader_arg_assignments_from_spec, **rader_arg_assignments_from_other_sources)
    
def _create_plink_links(bed_file_path, bim_file_path, fam_file_path, temp_dir):
    links_path_prefix = os.path.join(temp_dir,  'plink_%s_%s' % (os.path.basename(bed_file_path).split('.')[0], getuser()))
    safe_symlink(bed_file_path, links_path_prefix + '.bed')
    safe_symlink(bim_file_path, links_path_prefix + '.bim')
    _create_fam_link_or_fix_it(fam_file_path, links_path_prefix + '.fam')
    return links_path_prefix
    
def _create_fam_link_or_fix_it(original_fam_file_path, fixed_fam_file_path):
        
    '''
    A horrible hack, required to process UKBB's FAM files. This unfortunate situation is the result of the bad interface provided by the pandas-plink
    package, coupled with a format violation on UKBB's side. Apparently, UKBB's FAM file deviates from PLINK's FAM format by setting the 6th column
    (of costum values) to be non-numeric (in UKBB, they store the batch). It also appears that the pandas-plink module is not flexible enough to handle
    with such a deviation. As a result, we need to modify the 6th column of the FAM file before saving it in a temporary path, so that it could be
    processed by pandas-plink.
    '''
    
    if os.path.exists(fixed_fam_file_path):
        return

    fam_df = pd.read_csv(original_fam_file_path, sep = ' ', header = None)
    
    if np.issubdtype(fam_df.iloc[:, 5].dtype, np.number):
        # The FAM file seems ok. Can just create a link for it
        safe_symlink(original_fam_file_path, fixed_fam_file_path)
    else:
        
        try:
            fam_df.iloc[:, 5] = fam_df.iloc[:, 5].astype(float)
        except ValueError:
            # We are forced to override these values
            fam_df.iloc[:, 5] = np.nan
                    
        try:
            fam_df.to_csv(fixed_fam_file_path, sep = ' ', index = False, header = False)
        except OSError as e:
            if e.errno == 17:
                log('%s: already exists after all.' % fixed_fam_file_path)
            else:
                raise e
        
_GENOTYPING_FORMAT_TO_READER_TYPE = {
    'plink': PlinkGenotypingReader,
    'bgen': BgenGenotypingReader,
}
