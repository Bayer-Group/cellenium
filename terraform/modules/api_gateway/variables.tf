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


variable "vpc_id" {
  type = string
}

variable "nlb_name" {
  type    = string
  default = "cellenium-nlb"
}

variable "nlb_subnet_ids" {
  type    = list(string)
  default = []
}

variable "nlb_target_group_name" {
  type    = string
  default = "cellenium-api-nlb-tg"
}

variable "nlb_security_group_name" {
  type    = string
  default = "cellenium-api-nlb-sg"
}

variable "nlb_port" {
  type = number
}

variable "ec2_security_group_id" {
  type = string
}

variable "ec2_security_group_port" {
  type = number
}

variable "ec2_target_group_port" {
  type = number
}

variable "ec2_instance_id" {
  type = string
}

variable "api_name" {
  type    = string
  default = "cellenium-api"
}

variable "api_log_group_name" {
  type    = string
  default = "cellenium-api-log-group"
}

variable "api_authorizer_name" {
  type    = string
  default = "cellenium-api-authorizer"
}

variable "api_authorizer_audiences" {
  type    = list(string)
  default = []
}

variable "api_authorizer_issuer" {
  type = string
}

variable "api_domain_name" {
  type = string
}

variable "api_domain_name_certificate_arn" {
  type = string
}

variable "vpc_link_name" {
  type    = string
  default = "cellenium_api_vpc_link"
}

variable "vpc_link_subnet_ids" {
  type = list(string)
}

variable "vpc_link_security_group_name" {
  type    = string
  default = "cellenium_api_vpc_link_sg"
}

variable "route53_hosted_zone_id" {
  type = string
}