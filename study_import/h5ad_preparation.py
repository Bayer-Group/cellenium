from typing import List
import tqdm
import scanpy as sc
import pandas as pd
from anndata import AnnData
import logging

logging.basicConfig(level=logging.INFO)


def cellenium_settings(
        adata: AnnData,
        main_sample_attributes: List[str]
):
    d = adata.uns.get('cellenium', {})
    adata.uns['cellenium'] = d

    assert isinstance(main_sample_attributes, list)
    for a in main_sample_attributes:
        if a not in adata.obs.columns:
            raise Exception(f"main_sample_attributes: {a} not in observations dataframe")
    d['main_sample_attributes'] = main_sample_attributes


def calculate_differentially_expressed_genes(
        adata: AnnData,
        diffexp_attributes: List[str],
        layer: str,
        ngenes=100,
        diff_exp_min_group_expr=0.1,
        diff_exp_min_group_fc=0.5,
        diff_exp_max_notgroup_expr=1
):
    result_dataframes = []
    for diffexp_attribute in tqdm.tqdm(diffexp_attributes, desc='diff.exp. genes'):
        valid_attribute_group_check = (adata.obs[diffexp_attribute].value_counts() > 1)
        attr_values = valid_attribute_group_check.index[valid_attribute_group_check].tolist()

        sc.tl.rank_genes_groups(adata, diffexp_attribute, groups=attr_values, method='wilcoxon', layer=layer,
                                use_raw=False,
                                n_genes=ngenes)
        sc.tl.filter_rank_genes_groups(adata, min_in_group_fraction=diff_exp_min_group_expr,
                                       min_fold_change=diff_exp_min_group_fc,
                                       max_out_group_fraction=diff_exp_max_notgroup_expr, use_raw=False,
                                       key="rank_genes_groups", key_added="rank_genes_groups_filtered")
        for attr_value in attr_values:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups_filtered", group=attr_value)
            # remove filtered elements
            df = df[~df["names"].isnull()]
            df['ref_attr_value'] = attr_value
            df['cmp_attr_value'] = '_OTHERS_'
            df['attribute_name'] = diffexp_attribute
            result_dataframes.append(df)
    adata.uns.pop('rank_genes_groups', None)
    adata.uns.pop('rank_genes_groups_filtered', None)
    result_dataframe = pd.concat(result_dataframes, axis=0).reset_index(drop=True)
    logging.info("calculate_differentially_expressed_genes: found a list of genes for these attributes: %s",
                 result_dataframe['attribute_name'].unique().tolist())

    d = adata.uns.get('cellenium', {})
    adata.uns['cellenium'] = d
    d['differentially_expressed_genes'] = result_dataframe.copy()
    return result_dataframe
