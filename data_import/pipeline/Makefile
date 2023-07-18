.PHONY = build_lambda_layer build_lambda build_docker_image
.SECONDARY:

build_lambda_layer:
	set -e
	rm -rf lambda_build
	rm -rf lambda_layer.zip
	mkdir -p lambda_build/python
	docker run --platform linux/amd64 --volume $(shell pwd)/lambda:/opt/requirements --volume $(shell pwd)/lambda_build:/opt/lambda --rm python:3.10.12 /bin/bash -c "pip install -r /opt/requirements/lambda_requirements.txt -t /opt/lambda/python"
	cd lambda_build && zip -r ../lambda_layer.zip . && cd ..
	rm -rf lambda_build


build_lambda:
	set -e
	rm -rf lambda.zip
	cd lambda && zip -r ../lambda.zip . && cd ..

build_batch:
	set -e
	docker buildx build --platform linux/amd64 -t cellenium_data_import -f ./batch/Dockerfile ./..
