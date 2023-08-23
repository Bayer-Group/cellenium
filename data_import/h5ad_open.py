import io

import mudata
import scanpy as sc
from anndata import AnnData
from mudata import MuData
from smart_open import open


def h5ad_h5mu_read(filename) -> AnnData | MuData:
    open_param = filename
    if filename.startswith("s3:"):
        # AnnData's read_h5ad passes the "filename" parameter to h5py.File, which supports file-like objects in
        # addition to filename strings. It is able to read an AnnData file directly from S3 using the python
        # file-like object abstraction smart_open provides, however it seeks a lot and that causes read performance
        # to drop significantly. So we're copying the h5ad file into an in-memory file and read from there.
        s3_file_like_obj = open(filename, "rb")
        open_param = io.BytesIO(s3_file_like_obj.read())
        s3_file_like_obj.close()

    data = mudata.read_h5mu(open_param) if filename.endswith("h5mu") else sc.read_h5ad(open_param)

    if not isinstance(open_param, str):
        open_param.close()

    return data
