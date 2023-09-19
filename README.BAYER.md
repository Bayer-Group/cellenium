
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

```
sudo su -
...
```
