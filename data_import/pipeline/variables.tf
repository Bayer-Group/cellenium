variable "local_docker_host" {
  description = "The local docker host to use for building and pushing images, identify by using `docker context ls`"
  default = "unix:///var/run/docker.sock"
}

variable vpc_id {

}

variable terraform_remote_state_bucket {
  default = "cellenium-terraform-state"
}

variable "cellenium_batch_job_definition_name" {
  default = "cellenium-batch-job-definition"
}

variable "cellenium_batch_job_queue_name" {
  default = "cellenium-batch-job-queue"
}

variable "cellenium_batch_compute_environment_name" {
  default = "cellenium-batch-compute-environment"
}

variable "cellenium_batch_compute_environment_type" {
  default = "FARGATE"
}

variable "cellenium_db_secret_name" {
  default = "cellenium-db-secret"
}

variable "cellenium_db_user" {
  default   = "cellenium"
  sensitive = true
}

variable "cellenium_db_password" {
  default   = "cellenium"
  sensitive = true
}

variable "cellenium_db_ip" {
  sensitive = true
}

variable "cellenium_db_port" {
  sensitive = true
}

variable "cellenium_db_name" {
  sensitive = true
}

variable "cellenium_ec2_access_security_group_id" {
}

variable "cellenium_ec2_access_security_group_port" {
}