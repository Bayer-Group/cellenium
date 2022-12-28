# cellenium
a FAIR and scalable interactive visual analytics app for scRNA-Seq data

## Setting up
```bash
docker-compose up
docker exec -it cellenium_postgres_1 bash -c 'set -e; for f in /database_schema/*.sql; do echo "Processing $f"; psql --username postgres --host=localhost --echo-errors --set ON_ERROR_STOP=on --file=$f; done'
conda env create -f study_import/environment.yml
conda activate cellenium_import
(cd study_import && python masterdata.py)
# Run the public_data_run_scripts.ipynb notebook now.
(cd study_import && python cellenium_import.py)
```
