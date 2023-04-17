FROM python:3.10-slim

WORKDIR /usr/src/app

COPY study_import_requirements.txt ./
RUN pip install --no-cache-dir -r study_import_requirements.txt

COPY . .
