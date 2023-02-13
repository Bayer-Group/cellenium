.PHONY = reset_database normal_studydata_import test_studydata_import huge_studydata_import
.SECONDARY:

reset_database:
	docker-compose exec postgres bash -c 'set -e; for f in /database_schema/*.sql; do echo "Processing $$f"; psql --username postgres --host=localhost --echo-errors --set ON_ERROR_STOP=on --file=$$f; done'
	PYTHONPATH=$$(pwd)/data_import python data_import/masterdata.py
	rm -f scratch/*.h5ad.imported

scratch/%.h5ad: data_import/public_data/%.ipynb
	@echo jupyter notebook $< is expected to produce h5ad file $@ and $@.html
	PYTHONPATH=$$(pwd)/data_import jupyter nbconvert --execute --to html --stdout  $<  > $@.html

scratch/%.h5ad.imported: scratch/%.h5ad
	PYTHONPATH=$$(pwd)/data_import python data_import/study_import.py $< --analyze-database
	echo "done" > $@

test_studydata_import: scratch/pancreas_atlas_subset.h5ad.imported

normal_studydata_import: scratch/pancreas_atlas.h5ad.imported \
  scratch/blood_covid.h5ad.imported \
  scratch/lung_asthma.h5ad.imported \
  scratch/tabula_muris_senis_heart.h5ad.imported \
  scratch/tabula_muris_senis_liver.h5ad.imported \
  scratch/tabula_sapiens_kidney.h5ad.imported

huge_studydata_880kcells_33kgenes_import: scratch/heart_failure_reichart2022.h5ad.imported

huge_studydata_880kcells_50genes_import: scratch/heart_failure_reichart2022_gene_subset.h5ad.imported
