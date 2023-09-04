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