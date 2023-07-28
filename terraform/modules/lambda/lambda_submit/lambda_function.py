import json
import pathlib
import boto3
from urllib import parse
import os
import logging

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
    if not all(map(lambda var: os.environ.get(var, False), REQUIRED_ENV_VARIABLES)):
        raise RuntimeError(f"Missing required environment variable")

    bucket = event["Records"][0]["s3"]["bucket"]["name"]
    key = parse.unquote_plus(
        event["Records"][0]["s3"]["object"]["key"], encoding="utf-8"
    )

    logging.info(f"Processing {key}")

    engine = get_aws_db_engine()
    with engine.begin() as connection:
        rs = connection.execute(
            text(f"SELECT study_id FROM public.study WHERE import_file=:key"), dict(key=key)
        )
        study_ids = [r[0] for r in rs]
        logging.info(f"Study ids: {study_ids}")

    if pathlib.Path(key).suffix.lower() not in (".h5ad", ".h5mu") or len(study_ids) != 1:
        s3 = boto3.client("s3")
        logging.warning(f"Deleting {key}")
        s3.delete_object(Bucket=bucket, Key=key)
        return

    logging.info(f"Submitting job for {key}")
    batch_client = boto3.client("batch")
    response = batch_client.submit_job(
        jobDefinition=os.environ.get("BATCH_JOB_DEFINITION"),
        jobName="".join(e for e in key if e.isalnum()),
        jobQueue=os.environ.get("BATCH_JOB_QUEUE"),
        parameters={"filename": f"s3://{bucket}/{key}", "analyze-database": "--analyze-database"},
        tags={
            'studyId': str(study_ids[0])
        },
    )
    logging.info(f"Submitted job {response['jobId']} for {key}")
