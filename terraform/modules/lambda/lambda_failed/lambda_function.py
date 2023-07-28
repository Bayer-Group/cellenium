import json
import logging
import os

import boto3
from furl import furl
from sqlalchemy import create_engine
from sqlalchemy.sql import text

if len(logging.getLogger().handlers) > 0:
    # The Lambda environment pre-configures a handler logging to stderr. If a handler is already configured,
    # `.basicConfig` does not execute. Thus we set the level directly.
    logging.getLogger().setLevel(logging.INFO)
else:
    logging.basicConfig(level=logging.INFO)

REQUIRED_ENV_VARIABLES = ["AWS_DB_SECRET"]


def get_aws_db_engine():
    client = boto3.client(
        service_name="secretsmanager",
        region_name=os.environ.get("AWS_REGION", "eu-central-1"),
    )

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

    return create_engine(
        url=f"postgresql://{db_user}:{db_password}@{db_host}:{db_port}/{db}"
    )


def lambda_handler(event, context):
    if not all(map(lambda var: os.environ.get(var, False), REQUIRED_ENV_VARIABLES)):
        raise RuntimeError(f"Missing required environment variable")

    logging.info(event)
    logging.info(context)

    engine = get_aws_db_engine()

    detail = event.get("detail", dict())
    if detail.get("status") != "FAILED":
        return

    with engine.begin() as connection:

        study_id = detail.get("tags", dict()).get("study_id")

        file_url = None
        if detail.get("parameters", dict()).get("filename") is not None:
            file_url = detail.get("parameters", dict()).get("filename")
        elif detail.get("container", dict()).get("command") is not None:
            file_url = detail.get("container", dict()).get("command")[1]

        bucket = None
        key = None
        if file_url is not None:
            bucket = furl(file_url).host
            key = "/".join(furl(file_url).path.segments)

        if study_id is None:
            if key is None:
                return

            rs = connection.execute(
                text(f"SELECT study_id FROM public.study WHERE import_file=:file_name"), dict(file_name=key)
            )
            study_ids = [r[0] for r in rs]
            if len(study_ids) != 1:
                return
            study_id = study_ids[0]

        study_id = int(study_id)
        logging.info(f"Setting study {study_id} to failed")

        connection.execute(
            text("UPDATE public.study SET import_failed=true WHERE study_id=:study_id"), dict(study_id=study_id)
        )


    logs_client = boto3.client("logs")
    log_stream_name = detail.get("container", dict()).get("logStreamName")
    if log_stream_name is not None:
        try:
            events = logs_client.get_log_events(logStreamName=log_stream_name, logGroupName="/aws/batch/job")
            log = "\n".join([event["message"] for event in events.get("events")])

            next_token = events.get("nextForwardToken")
            while next_token is not None:
                events = logs_client.get_log_events(logStreamName=log_stream_name, logGroupName="/aws/batch/job")

                if events.get("nextForwardToken") == next_token:
                    next_token = None
                    continue
                else:
                    next_token = events.get("nextForwardToken")
                log += "\n" + "\n".join([event["message"] for event in events.get("events")])

            with engine.begin() as connection:
                connection.execute(
                    text("UPDATE public.study SET import_log=:log WHERE study_id=:study_id"),
                    dict(study_id=study_id, log=log)
                )

        except Exception as e:
            logging.error(e)

    if key is not None and bucket is not None:
        s3 = boto3.client("s3")
        logging.info(f"Deleting {key}")
        s3.delete_object(Bucket=bucket, Key=key)

    return "function ran", event