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

variable "stage" {
  type    = string
  default = "prod"
}


variable "docker_host" {
  type        = string
  description = "The local docker host to use for building and pushing images, identify by using `docker context ls`"
  default     = "unix:///var/run/docker.sock"
}

variable "vpc_id" {
  type = string
}

variable "subnet_ids" {
  type = list(string)
}

variable "ecr_repository_name" {
  type    = string
  default = "cellenium-study-import"
}

variable "job_definition_name" {
  type    = string
  default = "cellenium-study-import-jd"
}

variable "job_queue_name" {
  type    = string
  default = "cellenium-study-import-queue"
}

variable "compute_environment_name" {
  type    = string
  default = "cellenium-study-import-compute-environment"
}

variable "compute_environment_type" {
  type    = string
  default = "FARGATE"
}

variable "security_group_name" {
  type    = string
  default = "cellenium_study_import_batch_security_group"
}

variable "ec2_security_group_id" {
  type = string
}

variable "ec2_security_group_port" {
  type = string
}

variable "failed_batch_cloudwatch_rule_name" {
  type    = string
  default = "cellenium-study-import-failed-batch-rule"
}

variable "batch_import_role_name" {
  type    = string
  default = "cellenium-study-import-batch-role"
}

variable "manage_ecs_policy_name" {
  type    = string
  default = "cellenium_manage_ecs"
}


variable "batch_import_policy_name" {
  type    = string
  default = "cellenium-study-import-batch-policy"
}


variable "db_secret_id" {
  type = string
}
variable "db_secret_arn" {
  type = string
}

variable "bucket_arn" {
  type = string
}

variable "failed_lambda_arn" {
  type = string
}