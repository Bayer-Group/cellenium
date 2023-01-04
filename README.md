# cellenium
a FAIR and scalable interactive visual analytics app for scRNA-Seq data

## Setting up
```bash
docker-compose up
docker exec -it cellenium-postgres-1 bash -c 'set -e; for f in /database_schema/*.sql; do echo "Processing $f"; psql --username postgres --host=localhost --echo-errors --set ON_ERROR_STOP=on --file=$f; done'
conda env create -f data_import/environment.yml
conda activate cellenium_import
(cd data_import && python masterdata.py)
(cd data_import && python study_preparation.py)
(cd data_import && python study_import.py)
```
