import json
import logging
import os
import pathlib
from urllib import parse

import boto3
from sqlalchemy import create_engine
from sqlalchemy.sql import text

if len(logging.getLogger().handlers) > 0:
    # The Lambda environment pre-configures a handler logging to stderr. If a handler is already configured,
    # `.basicConfig` does not execute. Thus we set the level directly.
    logging.getLogger().setLevel(logging.INFO)
else:
    logging.basicConfig(level=logging.INFO)

REQUIRED_ENV_VARIABLES = ["BATCH_JOB_DEFINITION", "BATCH_JOB_QUEUE", "AWS_DB_SECRET"]


def get_aws_db_engine():
    client = boto3.client(
        service_name="secretsmanager",
        region_name=os.environ.get("AWS_REGION", "eu-central-1"),
    )
    logging.info(f"Getting secret {os.environ.get('AWS_DB_SECRET')}")
    secret_values = json.loads(
        client.get_secret_value(SecretId=os.environ.get("AWS_DB_SECRET"))[
            "SecretString"
        ]
    )
    db = secret_values["DB"]
    db_host = secret_values["IP"]
    db_port = secret_values["PORT"]
    db_user = secret_values["USER"]
    db_password = secret_values["PASSWORD"]

    logging.info(f"Connecting to {db_host}:{db_port}/{db}")

    return create_engine(
        url=f"postgresql://{db_user}:{db_password}@{db_host}:{db_port}/{db}"
    )


def lambda_handler(event, context):
    if not all(os.environ.get(var, False) for var in REQUIRED_ENV_VARIABLES):
        raise RuntimeError("Missing required environment variable")

    bucket = event["Records"][0]["s3"]["bucket"]["name"]
    key = parse.unquote_plus(
        event["Records"][0]["s3"]["object"]["key"], encoding="utf-8"
    )

    logging.info(f"Processing {key} {bucket}")
    if not key.startswith("s3://"):
        key = f"s3://{bucket}/{key}"

    engine = get_aws_db_engine()
    with engine.begin() as connection:
        rs = connection.execute(
            text("SELECT study_id, import_started, import_finished FROM public.study WHERE import_file=:key"), {"key": key}
        ).fetchall()
        study_ids = [r.study_id for r in rs]
        import_started = [r.import_started for r in rs]
        import_finished = [r.import_finished for r in rs]
        if len(study_ids) != 1:
            logging.error(f"Skipping S3 object {key}, filename mismatch")
            return
        logging.info(f"Study ID: {study_ids[0]}")
        if import_started[0] and not import_finished[0]:
            logging.error(f"Skipping S3 object {key}, it is currently being import")
            return

    if pathlib.Path(key).suffix.lower() not in (".h5ad", ".h5mu") or len(study_ids) != 1:
        s3 = boto3.client("s3")
        logging.warning(f"Deleting {key}")
        s3.delete_object(Bucket=bucket, Key=key)
        return

    job_name = "".join(e for e in key if e.isalnum())[:80]
    logging.info(f"Submitting job '{job_name}' for {key}")
    batch_client = boto3.client("batch")
    response = batch_client.submit_job(
        jobDefinition=os.environ.get("BATCH_JOB_DEFINITION"),
        jobName=job_name,
        jobQueue=os.environ.get("BATCH_JOB_QUEUE"),
        parameters={"filename": key},
        tags={
            'studyId': str(study_ids[0])
        },
    )
    logging.info(f"Submitted job {response['jobId']} for {key}")
