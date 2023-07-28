terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 4.16"
    }
    docker = {
      source  = "kreuzwerker/docker"
      version = "~>3.0.2"
    }
  }
}

provider "aws" {
  region = var.aws_region
}


data "aws_vpc" "vpc" {
  id = var.vpc_id
}

data "aws_caller_identity" "current" {}