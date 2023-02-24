from smart_open import open
import scanpy as sc
import io


def h5ad_read(filename):
    if filename.startswith('s3:'):
        # AnnData's read_h5ad passes the "filename" parameter to h5py.File, which supports file-like objects in
        # addition to filename strings. It is able to read an AnnData file directly from S3 using the python
        # file-like object abstraction smart_open provides, however it seeks a lot and that causes read performance
        # to drop significantly. So we're copying the h5ad file into an in-memory file and read from there.
        s3_file_like_obj = open(filename, 'rb')
        memory_file_like_obj = io.BytesIO(s3_file_like_obj.read())
        s3_file_like_obj.close()
        adata = sc.read_h5ad(memory_file_like_obj)
        memory_file_like_obj.close()
        return adata
    else:
        return sc.read_h5ad(filename)
