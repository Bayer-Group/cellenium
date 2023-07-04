FROM python:3.10-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends build-essential cmake

COPY ./pipeline/batch/requirements.txt .
RUN pip install --target=/opt/requirements -r ./requirements.txt


FROM python:3.10-slim AS runtime
COPY --from=build /opt/requirements/  /opt/requirements/
ENV PYTHONPATH="${PYTHONPATH}:/opt/requirements/"

WORKDIR /app
COPY study_import.py .
COPY h5ad_open.py .
COPY postgres_utils.py .

ENTRYPOINT ["python", "study_import.py"]

