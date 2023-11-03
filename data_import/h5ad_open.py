import shutil
import tempfile

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
        # to drop significantly. So we could copy the h5ad file into an in-memory file and read from there.
        # However, the h5mu library will take real file handles only. So lets resort to copying to a temp file
        # and use that.
        s3_file_like_obj = open(filename, "rb")
        fp = tempfile.NamedTemporaryFile(delete=False)
        shutil.copyfileobj(s3_file_like_obj, fp, 16 * 1024**2)
        s3_file_like_obj.close()
        fp.close()
        open_param = fp.name

    data = mudata.read_h5mu(open_param) if filename.endswith("h5mu") else sc.read_h5ad(open_param)

    if not isinstance(open_param, str):
        open_param.close()

    return data
