import scanpy as sc
import numpy as np
import h5ad_preparation as prep


def add_cellenium_settings(adata):
    prep.cellenium_settings(
        adata,
        main_sample_attributes=['cell_type', 'cell_states', 'disease', 'tissue']
    )


def download_prepare_h5ad():
    # Study with 881k cells and 1357m expression values
    # URL from https://cellxgene.cziscience.com/collections/e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb
    # TODO use sfaira for download (URL will expire) and align with general data handling code - this is just for testing.
    url = "https://corpora-data-prod.s3.amazonaws.com/aca04c31-6509-4409-8d96-5b0e5708948c/local.h5ad?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=ASIATLYQ5N5X6QH3GTMP%2F20230103%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20230103T091923Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEIj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJHMEUCIDu1a%2B5p5L%2FYYSctvKeAp9uPlSNP3QAN2UwbP%2Bck5XakAiEAwx2v8GXZwTGdfxtn1YG3PMFshkuZv7A5TxyPc1Fcopsq9AMI0f%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARABGgwyMzE0MjY4NDY1NzUiDArtz4L8MSBrqrW1WirIA7cGultk4R7FFTYzfJCLnJxxS%2FaA6wAzpF%2BfKn9E0fmWslV6UXztUCJEx5QXg%2FCIqHgEQ9WRkBO448JJONHQbaePkIYAo46%2BujqnG%2F1A%2F%2FFiZDMVEHyr%2FUxq0NT87%2Br4SPEoFlI4wymigm7NWOv%2B%2Fy7vfHmg0LLHJ%2FDZBFHJ949pAeN9naNYKw5YgYQ9E1Yde8CSbShto1yOY3sV0HClxWKp%2BIYBFs1HyzeCviAva3hHCVXg3J0%2Fa%2Bn7%2FrK7zezJaBLsatWTpvnXa1zmNeA6B1MIZK7tq32Blb%2FKnQevlozdpZJX9ANkI%2BGtE0LqsLzP3EyXaVUMdbqqdWD6%2B3k5cyA%2BTGVw2Loy23va9sH9ll%2Bw%2Berdp3P0iY%2FWWNaa0sye0SgDdRSVc2dRL6WOOusj6EQUJy6JIvJiPt7qLOGKS5oBa7HcIcv%2FAJwtlDYUrWCpt2EX7X61jmrlvjP562pVCI%2Bz9nrfHpb%2B3PCV3PqgfoYDrRg6GkEjjmOSWWtuFt9sDmB%2BLMZAJhFRx3kO3dcqDoBjmIfT1h%2FB10h68ceQejjtbz4TnIySmCfyUXiGo3pwJiKqAUKUoCUKiBrBZUWYnUaj9NDuR3xE%2FTDFs8%2BdBjqlAfE3IAV8sIXhwOSOb7FHJ5E1%2BalaXnHvrdkOfq3YbHjQeWRcJMQK6Z34fG4QFP6P7u%2Bb6HWxjo0M8ohlLF22ULnnSKmc2LM0Hlg0ACubhM5gUysTwGMimkHFBmHhtNpMNswaE1byGB4s1EoNlMjQu5zSwsiIa1425UFBKarx0KaDLP%2B6pCJ%2BMqiwpLJ6XrGrG0eRHJuXc%2BLUVUzPbml4RqwNWXl48g%3D%3D&X-Amz-Signature=0858c7f075fc24b68eba54702dc280c3dc9c0e94a31405a58c6fc2ae34e93f9d"
    adata = sc.read("../scratch/heart_failure_reichart2022_original.h5ad", backup_url=url)
    adata.layers['normalized'] = adata.X
    add_cellenium_settings(adata)
    #prep.calculate_differentially_expressed_genes(adata, ['celltype'], 'normalized')
    adata.write('../scratch/heart_failure_reichart2022.h5ad')

    query = np.array([s in ["MY0"] for s in adata.obs.cell_states])
    adata_subset = adata[query].copy()
    sc.pp.highly_variable_genes(
        adata_subset,
        n_top_genes=200,
        subset=True,
        layer='normalized'
    )
    adata_subset = adata_subset[query].copy()
    adata_subset = adata_subset[:, adata_subset.var_names].copy()
    adata_subset.layers['normalized'] = adata_subset.X
    add_cellenium_settings(adata_subset)
    prep.calculate_differentially_expressed_genes(adata_subset, ['disease'], 'normalized')
    adata_subset.write('../scratch/heart_failure_reichart2022_subset.h5ad')
