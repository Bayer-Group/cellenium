
### Deployment of webapp

`.env` file:

```
S3_BUCKET=850928066808-prod-cellenium-study-import
USER_STUDY_UPLOAD_CONFIGURED=true
AWS_REGION=eu-central-1
AUTH_HEADER_NAME=x-amzn-oidc-accesstoken
```

```
ubuntu@prod-c-cellenium2:/data/app/dtlab_cellenium_internal
git status
# we always use 'develop' branch while testing it here...
git pull
conda activate cellenium_import
conda env update -f data_import/environment.yml 
docker compose up -d --build
```

### Deployment of Lambda functions, API Gateway

`/data/app/dtlab_cellenium_internal/terraform/terraform.tfvars` file:

```
stage = "prod"

ec2_security_group_id                = "sg-0724ba548ae3e7ae1"
ec2_postgres_security_group_port     = 5001
ec2_postgraphile_security_group_port = 5000
vpc_id                               = "vpc-0ebca1a7e15b24110"
subnet_ids                           = ["subnet-0a5cab2045faac789"]

db_credentials = {
  username      = "postgres"
  password      = "...TODO SET ME..."
  ip            = "10.145.10.205"
  port          = 5001
  database_name = "postgres"
}

ec2_instance_id          = "i-05d1e0ae821a6716a"
route53_hosted_zone_id   = "Z04237652N01CSYQKIF4E"
nlb_subnet_ids           = ["subnet-0a5cab2045faac789"]
nlb_port                 = 5000
api_authorizer_issuer    = "https://sts.windows.net/fcb2b37b-5da0-466b-9b83-0014b67a7c78/"
api_authorizer_audiences = [
  "93f415d7-5aac-473e-b704-df25555b7177",
  "api://af54600e-373f-4993-a5ee-5a553da6fc7c",
  "af54600e-373f-4993-a5ee-5a553da6fc7c",
  "6da2e570-63dd-4882-af8b-6cfc1d0c94fc",
  "api://6da2e570-63dd-4882-af8b-6cfc1d0c94fc/.default",
  "api://6da2e570-63dd-4882-af8b-6cfc1d0c94fc",
  "1bebd04c-43e7-4a38-907c-7538252365b1",
]
api_domain_name                 = "cellllenium-api.chegenara.int.bayer.com"
api_domain_name_certificate_arn = "arn:aws:acm:eu-central-1:850928066808:certificate/e75b3732-7bb4-4258-9b07-4d26168e185e"
vpc_link_subnet_ids             = ["subnet-0a5cab2045faac789"]
```


```
sudo su -
cd /data/app/dtlab_cellenium_internal/terraform/
terraform apply
# to deploy a specific module only:
terraform apply --target=module.cellenium_study_import_lambda
```
