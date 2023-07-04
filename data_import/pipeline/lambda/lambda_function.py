import pathlib
import boto3
from urllib import parse
import os
import logging

logging.basicConfig(level=logging.WARNING)

client = boto3.client('batch')
s3 = boto3.client('s3')

REQUIRED_ENV_VARIABLES = ["BATCH_JOB_DEFINITION", "BATCH_JOB_QUEUE"]


def lambda_handler(event, context):
    if not all(map(lambda var: os.environ.get(var, False), REQUIRED_ENV_VARIABLES)):
        raise RuntimeError(f"Missing required environment variable")

    bucket = event['Records'][0]['s3']['bucket']['name']
    key = parse.unquote_plus(event['Records'][0]['s3']['object']['key'], encoding='utf-8')

    if pathlib.Path(key).suffix != '.h5ad':
        logging.warning(f"Deleting {key}")
        s3.delete_object(Bucket=bucket, Key=key)
        return

    response = client.submit_job(
        jobDefinition=os.environ.get("BATCH_JOB_DEFINITION"),
        jobName="".join(e for e in key if e.isalnum()),
        jobQueue=os.environ.get("BATCH_JOB_QUEUE"),
        parameters={
            "filename": f"s3://{bucket}/{key}",
            "analyze-database": "true"
        }
    )
    print(response)
