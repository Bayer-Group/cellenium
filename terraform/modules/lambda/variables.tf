variable "aws_region" {
  type    = string
  default = "eu-central-1"

  validation {
    condition = contains([
      "us-east-1",
      "us-east-2",
      "us-west-1",
      "us-west-2",
      "ca-central-1",
      "eu-west-1",
      "eu-central-1",
      "eu-west-2",
      "eu-west-3",
      "eu-north-1",
      "ap-northeast-1",
      "ap-northeast-2",
      "ap-southeast-1",
      "ap-southeast-2",
      "ap-south-1",
      "sa-east-1",
      "us-gov-west-1",
      "us-gov-east-1"
    ], var.aws_region)
    error_message = "Please provide a valid AWS region."
  }
}


variable "vpc_id" {
  type = string
}

variable "s3_bucket_name" {
  type = string
}

variable "s3_bucket_arn" {
  type = string
}

variable "s3_bucket_id" {
  type = string
}

variable "failed_batch_cloudwatch_rule_arn" {
  type = string
}

variable "ec2_security_group_id" {
  type = string
}

variable "ec2_security_group_port" {
  type = string
}

variable "batch_job_definition_name" {
  type = string
}

variable "batch_job_definition_arn" {
  type = string
}

variable "batch_queue_name" {
  type = string
}

variable "batch_queue_arn" {
  type = string
}

variable "db_secret_name" {
  type = string
}

variable "db_secret_arn" {
  type = string
}