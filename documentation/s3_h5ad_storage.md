# S3 for H5AD data storage (optional)

Cellenium uses h5ad files during study import (to postgres) and also interactively, when
calculating differentially expressed genes for a user-defined sample annotation group
and when finding genes whose expression is correlated to a query gene. Those operations require
the full data matrix, which is best obtained by reading the h5ad file.

In the default demo setup, h5ad files are stored in `./scratch` which is seen as `/h5ad_store`
in the postgres docker container.

As an optional feature, the `study_import.py` can be pointed to an S3 object instead of a
local filename. The data will be imported from S3, and the S3 path will be retained in
the postgres column `study.filname` for eventual download during the interactive use
cases described above.

All programs (`study_import.py`, the postgres docker container) will get their AWS credentials
either through environment variables or simply because they are run in an AWS environment
(e.g. EC2, ECS) which has the right IAM role set.

# User H5AD data upload

Further, cellenium can be configured to allow users to upload h5ad files through temporary
AWS credentials that allow to write in a prefix of the S3 bucket, using the normal AWS cli.
In addition to using S3 for h5ad storage, this involves setting up an IAM user (`cellenium-import`
in the example below), an IAM policy (`c-cellenium-bucket-input`) and an IAM role (`c-cellenium-importer`).

To give access to a user-scoped S3 prefix only, cellenium needs to get the current user's email
address. This is done with the JWT token that is also used for study-level permissions.

# AWS configuration

This combines S3 bucket and user data upload setup.

### create IAM user 'cellenium-import'

* IAM user name, e.g. `cellenium-import` in our example
* type: "programmatic access" only
* no permissions attached

### create IAM role for application

* e.g. `appserver` in our example
* assigned to the ECS container or EC2 instance that runs the docker-compose stack
* no permissions necessary (assigned in bucket policy below)

### create S3 bucket

* bucket name: e.g. `cellenium` in our example
* keep "Block all public access" setting
* enable service-side encryption as needed
* after creation, edit bucket:
* enable "Server access logging" as needed
* Permissions, paste Bucket Policy:

```
{
    "Version": "2012-10-17",
    "Id": "AllowAccess",
    "Statement": [
        {
            "Sid": "bucket_restrict",
            "Effect": "Allow",
            "Principal": {
                "AWS": [
                    "arn:aws:iam::123456789012:user/cellenium-import",
                    "arn:aws:iam::123456789012:role/appserver"
                ]
            },
            "Action": "s3:*",
            "Resource": [
                "arn:aws:s3:::cellenium",
                "arn:aws:s3:::cellenium/*"
            ]
        }
    ]
}
```

### create IAM policy and role for cellenium importer

This role can be used by using the generated temporary credentials (study admin page).
It grants access to the S3 bucket, `input/username/` prefix. The appserver can "assume" the role and
obtain temporary credentials.

steps:

* create policy `c-cellenium-bucket-input` (see json below)
* create role `c-cellenium-importer`, EC2 use case
    * attach policy `c-cellenium-bucket-input` to this role
    * description `Role that is assumed (STS) by users who can add input h5ad data in the cellenium bucket`
    * set Maximum session duration: 12 hours
    * in trust relationship tab, add "c-cellenium-importer-trust-relationship" statement (see json below)
    * the created role is used in the `S3_IMPORTER_ROLE_ARN` environment variable

"c-cellenium-bucket-input" policy:

```
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "ListObjectsInBucket",
            "Effect": "Allow",
            "Action": [
                "s3:ListBucket"
            ],
            "Resource": [
                "arn:aws:s3:::cellenium-test"
            ]
        },
        {
            "Sid": "AllObjectActions",
            "Effect": "Allow",
            "Action": "s3:*Object",
            "Resource": [
                "arn:aws:s3:::cellenium-test/input/*"
            ]
        }
    ]
}
```

"c-cellenium-importer-trust-relationship":

```
{
  "Effect": "Allow",
  "Principal": {
    "AWS": "arn:aws:iam::123456789012:role/appserver"
  },
  "Action": "sts:AssumeRole"
}

```

# Local development

You can obtain the role credentials for local development from e.g. an EC2 instance that has
the role (`appserver`) attached. See
https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html section "Retrieve security credentials
from instance metadata"

```
TOKEN=`curl -X PUT "http://169.254.169.254/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600"` && CREDS=`curl -H "X-aws-ec2-metadata-token: $TOKEN"  http://169.254.169.254/latest/meta-data/iam/security-credentials/c-appserver` && echo export AWS_ACCESS_KEY_ID=`jq -r  '.AccessKeyId' <<< "${CREDS}"` && echo export AWS_SECRET_ACCESS_KEY=\"`jq -r  '.SecretAccessKey' <<< "${CREDS}"`\" && echo export AWS_SESSION_TOKEN=\"`jq -r  '.Token' <<< "${CREDS}"`\"
```



