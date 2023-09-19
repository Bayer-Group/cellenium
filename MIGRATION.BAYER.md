
## cellenium v1:

export jupyter notebook in old codebase, has SQL statement for finding study IDs to export


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
```


## cellenium v2:

```
make -f Makefile.BAYER download
rm scratch/bayer_studies/*.imported
make -f Makefile.BAYER -k import_all
```