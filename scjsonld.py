import scanpy as sc
from anndata import AnnData
from typing import Union, Optional, Any, Mapping, Callable, Tuple, Sequence, Type
import numpy as np
from scipy.sparse import spmatrix
from numpy.random import RandomState
from louvain.VertexPartition import MutableVertexPartition
#import json

def read_10x_mtx(
        path,
        var_names='gene_symbols',
        make_unique=True,
        cache=False,
        gex_only=True):
    data = sc.read_10x_mtx(path,
                           var_names,
                           make_unique,
                           cache,
                           gex_only)
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['source']={'path': path,
                                    'var_names': var_names,
                                    'workflow_step': data.uns['workflow']['step']}
    return data

     
def filter_cells(
            data: AnnData,
            min_counts: Optional[int] = None,
            min_genes:  Optional[int] = None,
            max_counts: Optional[int] = None,
            max_genes:  Optional[int] = None,
            inplace: bool = True,
            copy: bool = False,
        ):
    sc.pp.filter_cells(data,
                    min_counts,
                    min_genes,
                    max_counts,
                    max_genes,
                    inplace,
                    copy)
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['filter_cells']={'min_counts': min_counts,
                                          'min_genes': min_genes,
                                          'max_counts': max_counts,
                                          'max_genes': min_genes,
                                          'workflow_step': data.uns['workflow']['step']}



def filter_genes(
            data: AnnData,
            min_counts: Optional[int] = None,
            min_cells:  Optional[int] = None,
            max_counts: Optional[int] = None,
            max_cells:  Optional[int] = None,
            inplace: bool = True,
            copy: bool = False,
        ):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['filter_genes']={'min_counts': min_counts,
                                          'min_cells': min_cells,
                                          'max_counts': max_counts,
                                          'max_cells': max_cells,
                                          'workflow_step': data.uns['workflow']['step']}
    sc.pp.filter_genes(data,
                    min_counts,
                    min_cells,
                    max_counts,
                    max_cells,
                    inplace,
                    copy)

def filter_obs_genes(
        data: AnnData,
        max_genes,
        ):
    tdata = data[data.obs['n_genes'] < max_genes, :]
    print(tdata.uns['workflow'].keys())
    if 'workflow' not in tdata.uns:
        if 'workflow' not in data.uns:
            tdata.uns['workflow']={'step':0}
        else:
            tdata.uns['workflow']=data.uns['workflow']
    tdata.uns['workflow']['step']=tdata.uns['workflow']['step']+1
    tdata.uns['workflow']['filter_obs_genes']={'max_genes': max_genes,
                                               'workflow_step': tdata.uns['workflow']['step']}
    return tdata


def filter_obs_percent_mito(
        data: AnnData,
        max_percent_mito,
        ):
    tdata = data[data.obs['percent_mito'] < max_percent_mito, :]
    if 'workflow' not in tdata.uns:
        if 'workflow' not in data.uns:
            tdata.uns['workflow']={'step':0}
        else:
            tdata.uns['workflow']=data.uns['workflow']
    tdata.uns['workflow']['step']=tdata.uns['workflow']['step']+1           
    tdata.uns['workflow']['filter_obs_percent_mito']={'max_percent_mito': max_percent_mito,
                                                      'workflow_step': data.uns['workflow']['step']}
                                                      
    return tdata

def normalize_per_cell(
        data: AnnData,
        counts_per_cell_after=None,
        counts_per_cell=None,
        key_n_counts=None,
        copy=False,
        layers=[],
        use_rep=None,
        min_counts=1):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1       
    data.uns['workflow']['normalize_per_cell']={'counts_per_cell_after': counts_per_cell_after,
                                                'counts_per_cell':counts_per_cell,
                                                'key_n_counts': key_n_counts,
                                                'use_rep': use_rep,
                                                'min_counts': min_counts,
                                                'workflow_step': data.uns['workflow']['step']}
                                                
    sc.pp.normalize_per_cell(data,
                             counts_per_cell_after,
                             counts_per_cell,
                             key_n_counts,
                             copy,
                             layers,
                             use_rep,
                             min_counts)


def log1p(
            data: Union[AnnData, np.ndarray, spmatrix],
            copy: bool = False,
            chunked: bool = False,
            chunk_size: Optional[int] = None,
        ):
    
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['log1p']={'executed': True,
                                   'chunked': chunked,
                                   'chunk_size': chunk_size,
                                   'workflow_step': data.uns['workflow']['step']}
                                   
    sc.pp.log1p(data, copy, chunked, chunk_size)

    
def highly_variable_genes(
            data,
            min_disp=None, max_disp=None,
            min_mean=None, max_mean=None,
            n_top_genes=None,
            n_bins=20,
            flavor='seurat',
            subset=False,
            inplace=True,
        ):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['high_variable_genes']={'min_disp': min_disp,
                                                 'max_disp': max_disp,
                                                 'min_mean': min_mean,
                                                 'max_mean': max_mean,
                                                 'n_top_genes': n_top_genes,
                                                 'n_bins': n_bins,
                                                 'flavor': flavor,
                                                 'subset': subset,
                                                 'workflow_step': data.uns['workflow']['step']}
                                                 
    sc.pp.highly_variable_genes(data,
                                min_disp,
                                max_disp,
                                min_mean,
                                max_mean,
                                n_top_genes,
                                n_bins,
                                flavor,
                                subset,
                                inplace)
                                           

def filter_highly_variable_genes(data: AnnData):
    tdata = data[:, data.var['highly_variable']]

    if 'workflow' not in tdata.uns:
        if 'workflow' not in data.uns:
            tdata.uns['workflow']={'step':0}
        else:
            tdata.uns['workflow']=data.uns['workflow']
    tdata.uns['workflow']['step']=tdata.uns['workflow']['step']+1            
    tdata.uns['workflow']['filter_highly_variable_genes']={'genes': tdata.var['highly_variable'],
                                                           'workflow_step': data.uns['workflow']['step']}

    return tdata


def regress_out(data, keys, n_jobs=None, copy=False):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['regress_out']={'keys': keys,
                                         'workflow_step': data.uns['workflow']['step']}
                                         
    sc.pp.regress_out(data, ['n_counts', 'percent_mito'], n_jobs, copy)

def scale(data, zero_center=True, max_value=None, copy=False):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['scale']={'zero_center': zero_center,
                                   'max_value': max_value,
                                   'workflow_step': data.uns['workflow']['step']}
                                   
    sc.pp.scale(data, zero_center, max_value, copy)

def pca(data: Union[AnnData, np.ndarray, spmatrix],
        n_comps: int = 50,
        zero_center: Optional[bool] = True,
        svd_solver: str = 'auto',
        random_state: int = 0,
        return_info: bool = False,
        use_highly_variable: Optional[bool] = None,
        dtype: str = 'float32',
        copy: bool = False,
        chunked: bool = False,
        chunk_size: Optional[int] = None):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['pca']={'n_comps': n_comps,
                                 'zero_center': zero_center,
                                 'svd_solver': svd_solver,
                                 'random_state': random_state,
                                 'return_info': return_info,
                                 'use_highly_variable': use_highly_variable,
                                 'dtype': dtype,
                                 'copy': copy,
                                 'chunked': chunked,
                                 'chunk_size': chunk_size,
                                 'workflow_step': data.uns['workflow']['step']}
                                 
    sc.tl.pca(data,
              n_comps,
              zero_center,
              svd_solver,
              random_state,
              return_info,
              use_highly_variable,
              dtype,
              copy,
              chunked,
              chunk_size)

def neighbors(data: AnnData,
              n_neighbors: int = 15,
              n_pcs: Optional[int] = None,
              use_rep: Optional[str] = None,
              knn: bool = True,
              random_state: Optional[Union[int, RandomState]] = 0,
              method: str = 'umap',
              metric: Union[str, Callable[[np.ndarray, np.ndarray], float]] = 'euclidean',
              metric_kwds: Mapping[str, Any] = {},
              copy: bool = False):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['neighbors']={'n_neighbors': n_neighbors,
                                       'n_pcs': n_pcs,
                                       'use_rep': use_rep,
                                       'knn': knn,
                                       'random_state': random_state,
                                       'method': method,
                                       'metric': metric,
                                       'metric_kwds': metric_kwds,
                                       'copy': copy,
                                       'workflow_step': data.uns['workflow']['step']}
                                       
    sc.pp.neighbors(data,
                    n_neighbors,
                    n_pcs,
                    use_rep,
                    knn,
                    random_state,
                    method,
                    metric,
                    metric_kwds,
                    copy)

def umap(data,
                min_dist=0.5,
                spread=1.0,
                n_components=2,
                maxiter=None,
                alpha=1.0,
                gamma=1.0,
                negative_sample_rate=5,
                init_pos='spectral',
                random_state=0,
                a=None,
                b=None,
                copy=False):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['umap']={'min_dist':min_dist,
                                  'spread': spread,
                                  'n_components': n_components,
                                  'maxiter': maxiter,
                                  'alpha': alpha,
                                  'gamma': gamma,
                                  'negative_sample_rate': negative_sample_rate,
                                  'init_pos': init_pos,
                                  'random_state': random_state,
                                  'a': a,
                                  'b': b,
                                  'copy': copy,
                                  'workflow_step': data.uns['workflow']['step']}
                                  
    sc.tl.umap(data, min_dist, spread, n_components, maxiter, alpha, gamma, negative_sample_rate, init_pos, random_state, a, b, copy)

def louvain(
            data: AnnData,
            resolution: Optional[float] = None,
            random_state: int = 0,
            restrict_to: Optional[Tuple[str, Sequence[str]]] = None,
            key_added: Optional[str] = None,
            adjacency: Optional[spmatrix] = None,
            flavor: str = 'vtraag',
            directed: bool = True,
            use_weights: bool = False,
            partition_type: Optional[Type[MutableVertexPartition]] = None,
            partition_kwargs: Optional[Mapping[str, Any]]=None,
            copy: bool = False,):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['louvain']={'resolution': resolution,
                                     'random_state': random_state,
                                     'restrict_to': restrict_to,
                                     'key_added': key_added,
                                     'flavor': flavor,
                                     'directed': directed,
                                     'use_weights': use_weights,
                                     'partition_type': partition_type,
                                     'partition_kwargs': partition_kwargs,
                                     'copy': copy,
                                     'workflow_step': data.uns['workflow']['step']}
                                     
    sc.tl.louvain(data,
                  resolution,
                  random_state,
                  restrict_to,
                  key_added,
                  adjacency,
                  flavor,
                  directed,
                  use_weights,
                  partition_type,
                  partition_kwargs,
                  copy)

def rank_genes_groups(
                data,
                groupby,
                use_raw=True,
                groups='all',
                reference='rest',
                n_genes=100,
                rankby_abs=False,
                key_added=None,
                copy=False,
                method='t-test_overestim_var',
                corr_method='benjamini-hochberg',
                **kwds):
    if 'workflow' not in data.uns:
        data.uns['workflow']={'step':0}
    data.uns['workflow']['step']=data.uns['workflow']['step']+1
    data.uns['workflow']['rank_genes_groups']={'groupby': groupby,
                                               'use_raw': use_raw,
                                               'groups': groups,
                                               'reference': reference,
                                               'n_genes': n_genes,
                                               'rankby_abs': rankby_abs,
                                               'key_added': key_added,
                                               'copy': copy,
                                               'method': method,
                                               'corr_method': corr_method,
                                               'workflow_step': data.uns['workflow']['step']}
                                               
    sc.tl.rank_genes_groups(data,
                            groupby,
                            use_raw,
                            groups,
                            reference,
                            n_genes,
                            rankby_abs,
                            key_added,
                            copy,
                            method,
                            corr_method,
                            **kwds)
    
#adata.write('./write/pbmc3k_withoutX.h5ad', compression='gzip')
