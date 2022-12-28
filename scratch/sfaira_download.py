import sfaira
import os
import scanpy as sc


basedir = '.'
datadir = os.path.join(basedir, 'raw')
metadir = os.path.join(basedir, 'meta')
cachedir = os.path.join(basedir, 'cache')
ds = sfaira.data.Universe(data_path=datadir, meta_path=metadir, cache_path=cachedir)

ds.subset(key="id", values=["homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x"]) #select only this dataset, we can get any sfaira dataset like that
ds.download()
ds.load(verbose=1)

ds.datasets['homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x'].streamline_metadata(schema="sfaira")  # convert the metadata annotation to the sfaira standard
adata = ds.datasets['homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x'].adata  # get the anndata object


