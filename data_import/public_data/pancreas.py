import scanpy as sc
from anndata import AnnData
import numpy as np
import h5ad_preparation as prep


def create_study_files():
    # pancreas url example
    url = "https://figshare.com/ndownloader/files/24539828"
    dat = prep.get_h5ad_from_url(url, "pancreas_atlas_original")

    # lung sfaira example
    #sfaira_id = "homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x"
    #adata = get_sfaira_h5ad(sfaira_id, basedir)

    # process the file to make sure expression, umap and metadata exist
    dat = prep.make_sparse(dat)
    dat = prep.make_norm_expression(dat, layer_from="counts")
    dat.layers.pop("counts")
    dat = prep.add_umap(dat)
    dat = prep.add_cellenium_settings(dat, ["celltype","tech"])
    dat.write(prep.basedir.joinpath("pancreas_atlas.h5ad"))

    # for testing, subsample the dataset & calculate differential expression
    query = np.array([s in ["celseq","smarter"] for s in dat.obs.tech])
    sc.pp.highly_variable_genes(dat, n_top_genes=200, batch_key="tech", subset=True, layer="norm_log_expression")
    dat_sub = dat[query].copy()
    dat_sub = dat_sub[:, dat.var_names].copy()

    # calculate differential expression
    dat_sub = prep.add_differential_expression_tables(dat_sub, ['celltype',"tech"], 'norm_log_expression')

    # write the final test file
    write_path = prep.basedir.joinpath("pancreas_atlas_subset.h5ad")
    dat_sub.write(write_path)
