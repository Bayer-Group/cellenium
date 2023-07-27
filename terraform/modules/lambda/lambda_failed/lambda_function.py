import json
import pathlib
import boto3
from urllib import parse
import os
import logging

from sqlalchemy import create_engine
from sqlalchemy.sql import text

logging.basicConfig(level=logging.WARNING)

client = boto3.client("batch")
s3 = boto3.client("s3")

REQUIRED_ENV_VARIABLES = ["AWS_DB_SECRET"]


def get_aws_db_engine():
    session = boto3.session.Session()
    client = session.client(
        service_name="secretsmanager",
        region_name=os.environ.get("AWS_REGION", "eu-central-1"),
    )

    secret_values = json.loads(
        client.get_secret_value(SecretId=os.environ.get("AWS_DB_SECRET"))[
            "SecretString"
        ]
    )
    db = secret_values["dbname"]
    db_host = secret_values["host"]
    db_port = secret_values["port"]
    db_user = secret_values["username"]
    db_password = secret_values["password"]

    return create_engine(
        url=f"postgresql://{db_user}:{db_password}@{db_host}:{db_port}/{db}"
    )


def lambda_handler(event, context):
    if not all(map(lambda var: os.environ.get(var, False), REQUIRED_ENV_VARIABLES)):
        raise RuntimeError(f"Missing required environment variable")


    # engine = get_aws_db_engine()
    # with engine.connect() as connection:
    #     rs = connection.execute(
    #         text("UPDATE public.study SET import_failed=true WHERE study_id=%s".format()), timeout=10
    #     )
    #     study_ids = [r[0] for r in rs]
    #
    # bucket = event["Records"][0]["s3"]["bucket"]["name"]
    # key = parse.unquote_plus(
    #     event["Records"][0]["s3"]["object"]["key"], encoding="utf-8"
    # )
    #
    # if pathlib.Path(key).suffix != ".h5ad":
    #     logging.warning(f"Deleting {key}")
    #     s3.delete_object(Bucket=bucket, Key=key)
    #     return
    #
    # response = client.submit_job(
    #     jobDefinition=os.environ.get("BATCH_JOB_DEFINITION"),
    #     jobName="".join(e for e in key if e.isalnum()),
    #     jobQueue=os.environ.get("BATCH_JOB_QUEUE"),
    #     parameters={"filename": f"s3://{bucket}/{key}", "analyze-database": "true"},
    # )
    return "function ran", event
