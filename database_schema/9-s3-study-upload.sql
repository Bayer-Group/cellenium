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



CREATE OR REPLACE FUNCTION create_study_upload(IN filetype text, IN study_name text)
    RETURNS json
    LANGUAGE plpython3u
    VOLATILE
AS
    $$
    import plpy
    import boto3
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

    if filetype not in [".h5ad", ".h5mu"]:
        raise Exception("filetype must be one of .h5ad or .h5mu")

    if 'S3_BUCKET' not in os.environ:
        raise Exception("S3_BUCKET environment variable not set")

    s3_client = boto3.client("s3")

    username = get_username()
    s3_prefix = f"/input/{username}/"
    study_name_cleaned = re.sub('[^0-9a-zA-Z]+', '*', study_name)
    s3_key = f"{s3_prefix}{study_name_cleaned}.{filetype}"

    # check if key_exists
    try:
        s3_client.head_object(Bucket=os.environ['S3_BUCKET'], Key=s3_key)
        s3_key = f"{s3_prefix}{study_name_cleaned}_{uuid.uuid4()}.{filetype}"
    except ClientError as e:
        if e.response['Error']['Code'] != '404':
            raise e

    response = s3_client.generate_presigned_post(
        os.environ['S3_BUCKET'],
        s3_key,
        Fields=None,
        Conditions=None,
        ExpiresIn=60**60*12, # 12 hours
    )
    plan = plpy.prepare("INSERT INTO study (study_name, filename, visible, import_started, import_file) VALUES ($1, $2, $3, $4, $5)", ["text", "text", "bool", "bool", "text"])
    plpy.execute(plan, [study_name, pathlib.Path(s3_key).name, False, False, s3_key])

    # The response contains the presigned URL and required fields
    return json.dumps(response)
    $$
