
## cellenium v1:

SQL statement for finding study IDs to export, and creating placeholder study records in cellenium 2:

```
select 'touch studies/' || id || '.id'
,'https://cellenium.chegenara.int.bayer.com/analysis/'|| id ,
               s.config ->> 'studyname',
       s.cell_count,
       s.config ->> 'user_groups',
       s.config ->> 'contact',
       'insert into study (study_id,study_name,import_file) values (' || id || ', ''migration placeholder'', ''' ||
       's3://850928066808-prod-cellenium-study-import/input/cellenium_v1_migration/' ||
       id || '__' || substr( replace(s.config ->> 'h5ad', 's3://cellenium-prod/input/', ''),
           position( '/' in replace(s.config ->> 'h5ad', 's3://cellenium-prod/input/', '') ) +1
           ) || ''');' insert_for_cellenium_v2
from study s
where (s.config ->> 'tax_id')::int in ( 9606,10090, 10116, 10117, 9541)
--and ( s.config ->> 'studytype' != 'scRNA')
and (s.config ->> 'studytype' is null or s.config ->> 'studytype' = 'scRNA')
and s.import_status='imported'
and (select count(1) from study_layer sl where sl.study_id = s.id) = 1
order by s.id;
```

Makefile

```
EXPORT_STUDIES := $(wildcard studies/*.id)
EXPORT_LOGS := $(EXPORT_STUDIES:.id=.html)

studies/%.html: studies/%.id
        echo working on study $$(basename $(patsubst %.id,%,$^) )
        ./run_notebook.sh $$(basename $(patsubst %.id,%,$^) )

export_all: $(EXPORT_LOGS)
```

run_notebook.sh

```
#!/bin/bash
echo running study $1
docker exec -it cellenium_celleniumworker_1 bash -c "STUDY_ID=$1 jupyter nbconvert --execute --to html --stdout notebooks/export_to_cellenium_v2.ipynb" >output_$1.html
#echo 'STUDY_ID=$1'
RET=$?
if [ "$RET" -eq "0" ]; then
        mv output_$1.html studies/$1.html;
else
        mv output_$1.html studies/$1.html_failed;
fi
```

```
cd /data/cellenium_legacy_export
make export_all

aws s3 cp s3://cellenium-prod/input_c2/8__wagner_fan_combined__final_UMI_counts_11Sep2020_layers_ensembl.h5ad s3://850928066808-prod-cellenium-study-import/input/cellenium_v1_migration/
```


## cellenium v2:

```
make -f Makefile.BAYER download
rm scratch/bayer_studies/*.imported
make -f Makefile.BAYER -k import_all
```