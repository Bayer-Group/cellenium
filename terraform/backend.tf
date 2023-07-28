terraform {
  backend "s3" {
    bucket = "cellenium-terraform-state"
    region = "eu-central-1"
  }
}