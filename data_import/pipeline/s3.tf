

resource "aws_s3_bucket" "cellenium_s3_bucket" {
  bucket = "${data.aws_caller_identity.current.account_id}-cellenium-study-upload"

  tags = {
    Name = "Cellenium Study Import"
  }

  lifecycle {
    prevent_destroy = true
  }
}