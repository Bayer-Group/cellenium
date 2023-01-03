import scanpy as sc
from anndata import AnnData
import h5ad_preparation as prep


basedir = "../../scratch"

# pancreas url example
url = "https://figshare.com/ndownloader/files/24539828"
dat = prep.get_h5ad_from_url(url, basedir, "pancreas_atlas")

# lung sfaira example
#sfaira_id = "homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x"
#adata = get_sfaira_h5ad(sfaira_id, basedir)

# process the file to make sure expression, umap and metadata exist
dat = prep.make_sparse(dat)
dat = prep.make_norm_expression(dat)
dat = prep.add_umap(dat)
dat = prep.add_cellenium_settings(dat, ["celltype","tech"])

# for testing, subsample the dataset & calculate differential expression
query = np.array([s in ["celseq","smarter"] for s in dat.obs.tech])
sc.pp.highly_variable_genes(dat, n_top_genes=200, batch_key="tech", subset=True)
dat_sub = dat[query].copy()
dat_sub = dat_sub[:, dat.var_names].copy()

# calculate differential expression
dat_sub = prep.add_differential_expression_tables(dat_sub, ['celltype',"tech"], 'norm_log_expression')

# write the final test file
write_path = basedir + "/pancreas_subset.h5ad"
dat_sub.write(write_path)