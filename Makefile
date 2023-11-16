src = data_import
black = black --line-length 140 --exclude '.*\.ipynb' $(src)
ruff = ruff ./$(pkg_src) --line-length 140 --select E,W,F,N,I,C,B,UP,PT,SIM,RUF --ignore E501,C901,B008,N815,N802,N803,SIM105

.PHONY = reset_database update_reference_data normal_studydata_import test_studydata_import huge_studydata_import atac_studydata_import cite_studydata_import format
.SECONDARY:

reset_database:
	docker compose exec postgres bash -c 'set -e; for f in /database_schema/*.sql; do echo "Processing $$f"; psql --username postgres --host=localhost --echo-errors --set ON_ERROR_STOP=on --file=$$f; done'
	PYTHONPATH=$$(pwd)/data_import python data_import/masterdata.py
	rm -f scratch/*.imported

update_reference_data:
	PYTHONPATH=$$(pwd)/data_import python data_import/masterdata.py

scratch/%: data_import/public_data/%.ipynb data_import/h5ad_preparation.py
	@echo jupyter notebook $< is expected to produce h5ad file $@ and $@.html
	PYTHONPATH=$$(pwd)/data_import jupyter nbconvert --execute --to html --stdout  $<  > $@.html

scratch/%.imported: scratch/%
	PYTHONPATH=$$(pwd)/data_import NUM_PROCESSES=1 python data_import/cellenium_cli.py study-import $< --analyze-database
	echo "done" > $@


test_studydata_import: scratch/pancreas_atlas_subset.h5ad.imported \
  scratch/brain3k_processed_subset.h5mu.imported

normal_studydata_import: scratch/pancreas_atlas.h5ad.imported \
  scratch/blood_covid.h5ad.imported \
  scratch/lung_asthma.h5ad.imported \
  scratch/tabula_muris_senis_heart.h5ad.imported \
  scratch/tabula_muris_senis_liver.h5ad.imported \
  scratch/tabula_sapiens_kidney.h5ad.imported

atac_studydata_import: scratch/brain3k_processed.h5mu.imported

cite_studydata_import: scratch/pbmc3k_processed.h5mu.imported

huge_studydata_880kcells_33kgenes_import: scratch/heart_failure_reichart2022.h5ad.imported

huge_studydata_880kcells_50genes_import: scratch/heart_failure_reichart2022_gene_subset.h5ad.imported

format:
	$(ruff) --fix
	$(black)

check-format:
	$(black) --check

lint:
	$(ruff) --format=github