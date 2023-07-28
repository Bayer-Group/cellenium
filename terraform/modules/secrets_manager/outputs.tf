output "secret_arn" {
  value = aws_secretsmanager_secret.db_secret.arn
}

output "secret_id" {
  value = aws_secretsmanager_secret.db_secret.id
}