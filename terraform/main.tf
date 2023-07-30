terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 4.16"
    }
  }
  required_version = ">= 1.5.0"
}

provider "aws" {
  region = "eu-central-1"
}


module "cellenium_study_import_s3" {
  source = "./modules/s3"
}


module "cellenium_study_import_lambda" {
  source = "./modules/lambda"

  vpc_id                           = var.vpc_id
  subnet_ids                       = var.subnet_ids
  batch_job_definition_name        = module.cellenium_study_import_batch.batch_job_definition_name
  batch_queue_name                 = module.cellenium_study_import_batch.batch_queue_name
  db_secret_id                     = module.cellenium_db_secret.secret_id
  db_secret_arn                    = module.cellenium_db_secret.secret_arn
  ec2_security_group_id            = var.ec2_security_group_id
  ec2_security_group_port          = var.ec2_security_group_port
  failed_batch_cloudwatch_rule_arn = module.cellenium_study_import_batch.failed_job_cloudwatch_rule_arn
  s3_bucket_arn                    = module.cellenium_study_import_s3.s3_bucket_arn
  s3_bucket_id                     = module.cellenium_study_import_s3.s3_bucket_id
  s3_bucket_name                   = module.cellenium_study_import_s3.s3_bucket_name
  batch_job_definition_arn         = module.cellenium_study_import_batch.batch_job_definition_arn
  batch_queue_arn                  = module.cellenium_study_import_batch.batch_queue_arn
}


module "cellenium_study_import_batch" {
  source = "./modules/batch"

  bucket_arn              = module.cellenium_study_import_s3.s3_bucket_arn
  db_secret_arn           = module.cellenium_db_secret.secret_arn
  db_secret_id            = module.cellenium_db_secret.secret_id
  ec2_security_group_id   = var.ec2_security_group_id
  ec2_security_group_port = var.ec2_security_group_port
  failed_lambda_arn       = module.cellenium_study_import_lambda.failed_lambda_arn
  vpc_id                  = var.vpc_id
  subnet_ids              = var.subnet_ids
  docker_host             = var.docker_host
}


module "cellenium_db_secret" {
  source = "./modules/secrets_manager"

  db_name  = var.db_credentials.database_name
  ip       = var.db_credentials.ip
  password = var.db_credentials.password
  port     = var.db_credentials.port
  user     = var.db_credentials.username
}