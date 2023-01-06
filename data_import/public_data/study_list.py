import scanpy as sc
from anndata import AnnData
import numpy as np
import pandas as pd
import h5ad_preparation as prep
from pathlib import Path


def create_study_files(sep = "\t"):
    
    study_list_filename = "public_data_links.txt"
    
    filedir = Path().parent.parent.joinpath('public_data').resolve()
    localfile = file.joinpath(f"{study_list_filename}")
    ll = pd.read_csv(path_to_study_table, sep = sep)
    
    for index, row in ll.head().iterrows():
        dat = prep.get_h5ad_from_url(row["Data.URL"], row["Short.Name"])
        dat = prep.make_sparse(dat)
        dat = prep.make_norm_expression(dat, layer_from="counts")
        dat.layers.pop("counts")
        dat = prep.add_umap(dat)
        dat = prep.add_cellenium_settings(dat, row["title"],
                                     row["description"], 
                                     row["Taxonomy.ID"], 
                                     row["Tissue.NCIT"], 
                                     row["Disease.MESH.ID"], 
                                     row["Pubmed.ID"], 
                                     row["Data.URL"])
        
        # TODO: subsetting, detecting cell type attributes, differential expression on automatically determined attributes, obtaining a common ontology cell type by fuzzy matching
        
        # write the final test file
        write_path = prep.basedir.joinpath(row["Short.Name"] + "_cellenium.h5ad")
        dat.write(write_path)
