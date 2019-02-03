import scanpy as sc
from anndata import AnnData
from typing import Optional

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
        data.uns['workflow']={}
    data.uns['workflow']['source']={'path': path,
                                    'var_names': var_names}
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
        data.uns['workflow']={}
    data.uns['workflow']['filter_cells']={'min_counts': min_counts,
                                          'min_genes': min_genes,
                                          'max_counts': max_counts,
                                          'max_genes': min_genes}


def filter_genes(
            data: AnnData,
            min_counts: Optional[int] = None,
            min_cells:  Optional[int] = None,
            max_counts: Optional[int] = None,
            max_cells:  Optional[int] = None,
            inplace: bool = True,
            copy: bool = False,
        ):
    sc.pp.filter_genes(data,
                    min_counts,
                    min_cells,
                    max_counts,
                    max_cells,
                    inplace,
                    copy)
    if 'workflow' not in data.uns:
        data.uns['workflow']={}
    data.uns['workflow']['filter_genes']={'min_counts': min_counts,
                                          'min_cells': min_cells,
                                          'max_counts': max_counts,
                                          'max_cells': max_cells}

def filter_obs_genes(
        data: AnnData,
        max_genes,
        ):
    tdata = data[data.obs['n_genes'] < max_genes, :]
    if 'workflow' not in tdata.uns:
        tdata.uns['workflow']={}
    tdata.uns['workflow']['filter_obs_genes']={'max_genes': max_genes}
    return tdata


def filter_obs_percent_mito(
        data: AnnData,
        max_percent_mito,
        ):
    tdata = data[data.obs['percent_mito'] < max_percent_mito, :]
    if 'workflow' not in tdata.uns:
        tdata.uns['workflow']={}
    tdata.uns['workflow']['filter_obs_genes']={'max_percent_mito': max_percent_mito}
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
    sc.pp.normalize_per_cell(data,
                             counts_per_cell_after,
                             counts_per_cell,
                             key_n_counts,
                             copy,
                             layers,
                             use_rep,
                             min_counts)
    if 'workflow' not in data.uns:
        data.uns['workflow']={}
    data.uns['workflow']['normalize_per_cell']={'counts_per_cell_after': counts_per_cell_after,
                                                'counts_per_cell':counts_per_cell,
                                                'key_n_counts': key_n_counts,
                                                'use_rep': use_rep,
                                                'min_counts': min_counts}
