variable "stage" {
  type    = string
  default = "prod"
}

variable "docker_host" {
  type        = string
  description = "The local docker host to use for building and pushing images, identify by using `docker context ls`"
  default     = "unix:///var/run/docker.sock"
}

variable "ec2_security_group_id" {
  type = string
}

variable "ec2_security_group_port" {
  type = number
}

variable "ec2_instance_id" {
  type = string
}

variable "vpc_id" {
  type = string
}

variable "subnet_ids" {
  type = list(string)
}

variable "db_credentials" {
  type = map(string)

  validation {
    condition = alltrue(
      [for key in keys(var.db_credentials) : contains(["ip", "port", "database_name", "password", "username"], key)]
    )
    error_message = "db_credentials must contain the following keys: ip, port, database_name, password, username"
  }
}

variable "route53_hosted_zone_id" {
  type = string
}

variable "nlb_subnet_ids" {
  type    = list(string)
  default = []
}

variable "nlb_port" {
  type = number
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

variable "vpc_link_subnet_ids" {
  type = list(string)
}