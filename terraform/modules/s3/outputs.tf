

output "s3_bucket_name" {
  value = "${data.aws_caller_identity.current.account_id}-${var.bucket_name}"
}

output "s3_bucket_arn" {
  value = aws_s3_bucket.cellenium_s3_bucket.arn
}

output "s3_bucket_id" {
  value = aws_s3_bucket.cellenium_s3_bucket.id
}