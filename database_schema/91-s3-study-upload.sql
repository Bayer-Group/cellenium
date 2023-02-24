-- for user upload of h5ad study data - see docs/s3_h5ad_storage.md for configuration overview

CREATE OR REPLACE FUNCTION user_study_upload_configured()
    RETURNS boolean
    LANGUAGE plpython3u
    IMMUTABLE
AS
$$

import os

return True if os.getenv('S3_IMPORTER_ROLE_ARN') else False
$$


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

