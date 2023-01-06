# cellenium
a FAIR and scalable interactive visual analytics app for scRNA-Seq data

## Setting up
```bash
docker-compose up
conda env create -f data_import/environment.yml
conda activate cellenium_import
# 'test_studydata' should contain data to cover all application features, but is small enough to be imported in a few minutes
make reset_database test_studydata_import
# 'normal_studydata': real life studies (i.e. with full amount of cells and genes)
make normal_studydata_import
```
