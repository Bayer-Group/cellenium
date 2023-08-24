-- for user upload of h5ad study data - see docs/s3_h5ad_storage.md for configuration overview

CREATE OR REPLACE FUNCTION user_study_upload_configured()
    RETURNS boolean
    LANGUAGE plpython3u
    IMMUTABLE
AS
$$

import os
return True if os.getenv('USER_STUDY_UPLOAD_CONFIGURED') else False
$$;


CREATE OR REPLACE FUNCTION create_s3_temp_credentials()
    RETURNS text[]
    LANGUAGE plpython3u
    VOLATILE
AS
$$


import plpy
import boto3
import os
from jose import jwt


def get_token():
    r = plpy.execute("SELECT current_setting('postgraphile.auth_header_value', TRUE)::VARCHAR token")
    return [row for row in r][0]['token']


def get_username():
    claims = jwt.get_unverified_claims(get_token())
    return claims['unique_name']


def aws_credentials_for_bucket_user_prefix():
    # return ["AccessKeyId", "SecretAccessKey", "SessionToken", "username", "s3_prefix"]

    username = get_username()
    s3_prefix = f"{os.environ['S3_BUCKET']}/input/{username}/"
    sts_client = boto3.client('sts')
    assumed_role_object = sts_client.assume_role(
        RoleArn=os.environ['S3_IMPORTER_ROLE_ARN'],
        RoleSessionName=username,
        Policy="""{
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Sid": "ListObjectsInBucket",
                    "Effect": "Allow",
                    "Action": [
                        "s3:ListBucket"
                    ],
                    "Resource": [
                        "arn:aws:s3:::_S3_BUCKET_"
                    ]
                },
                {
                    "Sid": "AllObjectActions",
                    "Effect": "Allow",
                    "Action": "s3:*Object",
                    "Resource":
                        "arn:aws:s3:::_S3_PREFIX_*"

                }
            ]
        }""".replace('_S3_BUCKET_', os.environ['S3_BUCKET']).replace('_S3_PREFIX_', s3_prefix),
        DurationSeconds=12 * 60 * 60)
    return [
        assumed_role_object['Credentials']['AccessKeyId'],
        assumed_role_object['Credentials']['SecretAccessKey'],
        assumed_role_object['Credentials']['SessionToken'],
        username,
        s3_prefix
    ]


return aws_credentials_for_bucket_user_prefix()
$$;



CREATE OR REPLACE FUNCTION create_study_upload(IN filename text)
    RETURNS json
    LANGUAGE plpython3u
    VOLATILE
AS
    $$
    import plpy
    import boto3
    from botocore.client import Config
    from botocore.exceptions import ClientError
    from jose import jwt
    import pathlib
    import json
    import os
    import uuid
    import re

    def get_token():
        r = plpy.execute("SELECT current_setting('postgraphile.auth_header_value', TRUE)::VARCHAR token")
        return [row for row in r][0]['token']

    def get_username():
        try:
            claims = jwt.get_unverified_claims(get_token())
            return claims['unique_name']
        except Exception as e:
            return "anonymous"

    if not (filename.endswith(".h5ad") or filename.endswith(".h5mu")):
        raise Exception("filename must have .h5ad or .h5mu suffix")

    if 'S3_BUCKET' not in os.environ:
        raise Exception("S3_BUCKET environment variable not set")
    if 'AWS_REGION' not in os.environ:
        raise Exception("AWS_REGION environment variable not set")

    kwargs = dict()
    if all([key in os.environ for key in ["AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"]]):
        kwargs = dict(aws_access_key_id=os.environ['AWS_ACCESS_KEY_ID'],
                      aws_secret_access_key=os.environ['AWS_SECRET_ACCESS_KEY'],
                      aws_session_token=os.environ['AWS_SESSION_TOKEN'])

    s3_client = boto3.client("s3", config=Config(signature_version='s3v4'), region_name=os.environ['AWS_REGION'], **kwargs)

    username = get_username()
    cleaned_username = re.sub('[^0-9a-zA-Z]+', '', username)
    s3_prefix = f"input/{cleaned_username}/"
    filename_cleaned = re.sub('[^0-9a-zA-Z]+', '_', filename)
    if len(filename_cleaned) == 0:
        filename_cleaned = str(uuid.uuid4())
    s3_key = f"{s3_prefix}{filename_cleaned}"

    # check if key exists, and make it unique with a random string before the file suffix
    try:
        s3_client.head_object(Bucket=os.environ['S3_BUCKET'], Key=s3_key)
        s3_key = f"{s3_prefix}{filename_cleaned[:-5]}_{uuid.uuid4()}.{filename_cleaned[-4:]}"
    except ClientError as e:
        if e.response['Error']['Code'] != '404':
            raise e
    response = s3_client.generate_presigned_post(
        os.environ['S3_BUCKET'],
        s3_key,
        Fields=None,
        Conditions=None,
        ExpiresIn=60 * 60 * 12,  # 12 hours
    )

    if username == "anonymous":
        user_access = None
    else:
        user_access = [username]
    plan = plpy.prepare(
        "INSERT INTO study (study_name, filename, visible, import_started, import_file, reader_permissions, admin_permissions) VALUES ($1, $2, $3, $4, $5, $6, $7)",
        ["text", "text", "bool", "bool", "text", "text[]", "text[]"])
    plpy.execute(plan, ['(name will be obtained from file metadata)', pathlib.Path(s3_key).name, False, False, f"s3://{os.environ['S3_BUCKET']}/{s3_key}", user_access, user_access])
    # The response contains the presigned URL and required fields
    return json.dumps(response)
    $$;
grant insert, select, update, delete on study to postgraphile;
grant usage, select ON sequence study_study_id_seq TO postgraphile;