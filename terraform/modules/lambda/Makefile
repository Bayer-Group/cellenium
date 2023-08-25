.PHONY = build_lambda_layer build_lambda_submit build_lambda_failed
.SECONDARY:

build_lambda_layer:
	set -e
	rm -rf lambda_build
	rm -rf lambda_layer.zip
	mkdir -p lambda_build/python
	docker run --platform linux/amd64 --volume $(shell pwd)/lambda_layer:/opt/requirements --volume $(shell pwd)/lambda_build:/opt/lambda --rm python:3.10.12 /bin/bash -c "pip install -r /opt/requirements/requirements.txt -t /opt/lambda/python"
	cd lambda_build && zip -r ../lambda_layer.zip . && cd ..
	rm -rf lambda_build


build_lambda_submit:
	set -e
	rm -rf lambda.zip
	cd lambda_submit && zip -r ../lambda_submit.zip . && cd ..

build_lambda_failed:
	set -e
	rm -rf lambda.zip
	cd lambda_failed && zip -r ../lambda_failed.zip . && cd ..