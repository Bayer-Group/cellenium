terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 4.16"
    }
    docker = {
      source  = "kreuzwerker/docker"
      version = "3.0.2"
    }
  }

#  backend "s3" {
#    bucket = var.terraform_remote_state_bucket
#    region = "eu-central-1"
#  }

  required_version = ">= 1.2.0"
}


# configure docker provider
provider "docker" {
  host = var.local_docker_host
  registry_auth {
    address  = data.aws_ecr_authorization_token.token.proxy_endpoint
    username = data.aws_ecr_authorization_token.token.user_name
    password = data.aws_ecr_authorization_token.token.password
  }
}

provider "aws" {
  region  = "eu-central-1"
}

data "aws_vpc" "vpc" {
  id = var.vpc_id
}

data "aws_subnets" "default_vpc_subnets" {
  filter {
    name   = "vpc-id"
    values = [data.aws_vpc.vpc.id]
  }
}


data "aws_caller_identity" "current" {}

# get authorization credentials to push to ecr
data "aws_ecr_authorization_token" "token" {}


resource "aws_secretsmanager_secret" "db_secret" {
  name = var.cellenium_db_secret_name
  recovery_window_in_days = 0
}

resource "aws_secretsmanager_secret_version" "db_secret_version" {
  secret_id     = aws_secretsmanager_secret.db_secret.id
  secret_string = jsonencode(tomap({
    IP       = var.cellenium_db_ip
    PORT     = var.cellenium_db_port
    USER     = var.cellenium_db_user
    PASSWORD = var.cellenium_db_password
    DB       = var.cellenium_db_name
  }))
}