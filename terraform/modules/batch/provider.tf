terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 4.16"
    }
  }
}

provider "aws" {
  region = var.aws_region
}


# get authorization credentials to push to ecr
data "aws_ecr_authorization_token" "token" {}

data "aws_vpc" "vpc" {
  id = var.vpc_id
}

data "aws_caller_identity" "current" {}