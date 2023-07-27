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

variable "secret_name" {
  type = string
  default = "cellenium-db-credentials"
}

variable "user" {
  type = string
  sensitive = true
}

variable "password" {
  type = string
  sensitive = true
}

variable "ip" {
  type = string
  sensitive = true
}

variable "port" {
  type = number
  sensitive = true
}

variable "db_name" {
  type = string
  sensitive = true
}