

resource "aws_s3_bucket" "cellenium_s3_bucket" {
  bucket = "${data.aws_caller_identity.current.account_id}-${var.bucket_name}"

  lifecycle {
    prevent_destroy = true
  }

}