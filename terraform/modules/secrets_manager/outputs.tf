output "secret_arn" {
  value = aws_secretsmanager_secret.db_secret.arn
}

output "secret_name" {
  value = aws_secretsmanager_secret.db_secret.name
}