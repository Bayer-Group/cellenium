terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = ">= 4.67"
    }
    docker = {
      source  = "kreuzwerker/docker"
      version = "3.0.2"
    }
    null = {
      source  = "hashicorp/null"
      version = ">= 3.0"
    }
  }
  required_version = ">= 1.5.0"
}

provider "aws" {
  region = "eu-central-1"
}


module "cellenium_study_import_s3" {
  source = "./modules/s3"
  stage  = var.stage
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
  ec2_security_group_port          = var.ec2_postgres_security_group_port
  failed_batch_cloudwatch_rule_arn = module.cellenium_study_import_batch.failed_job_cloudwatch_rule_arn
  s3_bucket_arn                    = module.cellenium_study_import_s3.s3_bucket_arn
  s3_bucket_id                     = module.cellenium_study_import_s3.s3_bucket_id
  s3_bucket_name                   = module.cellenium_study_import_s3.s3_bucket_name
  batch_job_definition_arn         = module.cellenium_study_import_batch.batch_job_definition_arn
  batch_queue_arn                  = module.cellenium_study_import_batch.batch_queue_arn
  stage                            = var.stage
}


module "cellenium_study_import_batch" {
  source = "./modules/batch"

  bucket_arn              = module.cellenium_study_import_s3.s3_bucket_arn
  db_secret_arn           = module.cellenium_db_secret.secret_arn
  db_secret_id            = module.cellenium_db_secret.secret_id
  ec2_security_group_id   = var.ec2_security_group_id
  ec2_security_group_port = var.ec2_postgres_security_group_port
  failed_lambda_arn       = module.cellenium_study_import_lambda.failed_lambda_arn
  vpc_id                  = var.vpc_id
  subnet_ids              = var.subnet_ids
  docker_host             = var.docker_host
  stage                   = var.stage
}


module "cellenium_db_secret" {
  source = "./modules/secrets_manager"

  db_name  = var.db_credentials.database_name
  ip       = var.db_credentials.ip
  password = var.db_credentials.password
  port     = var.db_credentials.port
  user     = var.db_credentials.username
}


module "api_gateway" {
  source = "./modules/api_gateway"

  stage                  = var.stage
  vpc_id                 = var.vpc_id
  route53_hosted_zone_id = var.route53_hosted_zone_id

  api_authorizer_issuer           = var.api_authorizer_issuer
  api_authorizer_audiences        = var.api_authorizer_audiences
  api_domain_name                 = var.api_domain_name
  api_domain_name_certificate_arn = var.api_domain_name_certificate_arn
  ec2_security_group_id           = var.ec2_security_group_id
  ec2_security_group_port         = var.ec2_postgraphile_security_group_port
  ec2_target_group_port           = var.ec2_postgraphile_security_group_port
  vpc_link_subnet_ids             = var.vpc_link_subnet_ids

  ec2_instance_id = var.ec2_instance_id
  nlb_port        = var.nlb_port
  nlb_subnet_ids  = var.nlb_subnet_ids
}