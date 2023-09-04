terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = ">= 4.67"
    }
  }
  required_version = ">= 1.5.0"
}


# get authorization credentials to push to ecr
data "aws_ecr_authorization_token" "token" {}

data "aws_caller_identity" "current" {}