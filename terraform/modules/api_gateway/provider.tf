data "aws_vpc" "vpc" {
  id = var.vpc_id
}

data "aws_caller_identity" "current" {}