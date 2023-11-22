FROM --platform=linux/amd64 python:3.10 AS build
RUN rm -rf /var/lib/apt/lists/* && apt-get clean -y && apt-get update -y && apt-get upgrade -y && apt-get install -y --fix-missing build-essential cmake

COPY ./data_import/study_import_requirements.txt .
RUN pip install --target=/opt/requirements -r ./study_import_requirements.txt


FROM --platform=linux/amd64 python:3.10-slim AS runtime
COPY --from=build /opt/requirements/  /opt/requirements/
ENV PYTHONPATH="${PYTHONPATH}:/app:/opt/requirements/"

WORKDIR /app
COPY ./data_import/common.py .
COPY ./data_import/postgres_utils.py .
COPY ./data_import/h5ad_open.py .
COPY ./data_import/h5ad_preparation.py .
COPY ./data_import/study_import.py .
COPY ./data_import/huge_palette.py .
COPY ./data_import/differential_expression_calculation.py .
COPY ./data_import/cellenium_cli.py .
ADD ./data_import/Density_Sampling /app/Density_Sampling

ENTRYPOINT ["python", "cellenium_cli.py"]

