resource "aws_secretsmanager_secret" "db_secret" {
  name                    = var.secret_name
  recovery_window_in_days = 0
}

resource "aws_secretsmanager_secret_version" "db_secret_version" {
  secret_id     = aws_secretsmanager_secret.db_secret.id
  secret_string = jsonencode(tomap({
    IP       = var.ip
    PORT     = var.port
    USER     = var.user
    PASSWORD = var.password
    DB       = var.db_name
  }))
}