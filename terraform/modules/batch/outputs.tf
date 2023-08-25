

output "batch_job_definition_name" {
  value = aws_batch_job_definition.cellenium_study_import_job_definition.name
}

output "batch_job_definition_arn" {
  value = aws_batch_job_definition.cellenium_study_import_job_definition.arn
}

output "batch_queue_name" {
  value = aws_batch_job_queue.cellenium_study_import_job_queue.name
}

output "batch_queue_arn" {
  value = aws_batch_job_queue.cellenium_study_import_job_queue.arn
}

output "failed_job_cloudwatch_rule_arn" {
  value = aws_cloudwatch_event_rule.failed_job_cloudwatch_event_rule.arn
}