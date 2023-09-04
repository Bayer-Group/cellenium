resource "null_resource" "lambda_layer_zip" {
  provisioner "local-exec" {
    command     = "make build_lambda_layer"
    working_dir = path.module
  }

  triggers = {
    source_code_hash = filebase64sha256("${path.module}/lambda_layer/requirements.txt")
    zip_file         = fileexists("${path.module}/lambda_layer.zip") ? "1" : uuid()
  }
}

data "archive_file" "lambda_submit_zip" {
  type        = "zip"
  source_file = "${path.module}/lambda_submit/lambda_function.py"
  output_path = "${path.module}/lambda_submit.zip"
}

data "archive_file" "lambda_failed_zip" {
  type        = "zip"
  source_file = "${path.module}/lambda_failed/lambda_function.py"
  output_path = "${path.module}/lambda_failed.zip"
}

resource "aws_security_group" "cellenium_study_import_lambda_security_group" {
  name        = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.lambda_security_group}"
  description = "Security group for cellenium study import with AWS Batch"
  vpc_id      = var.vpc_id
}

resource "aws_security_group_rule" "cellenium_study_import_lambda_egress_security_group_rule" {
  type                     = "egress"
  from_port                = var.ec2_security_group_port
  to_port                  = var.ec2_security_group_port
  protocol                 = "tcp"
  source_security_group_id = var.ec2_security_group_id
  security_group_id        = aws_security_group.cellenium_study_import_lambda_security_group.id
}


resource "aws_security_group_rule" "cellenium_study_import_lambda_egress_internet_security_group_rule" {
  type              = "egress"
  from_port         = 0
  to_port           = 65535
  protocol          = "tcp"
  cidr_blocks       = ["0.0.0.0/0"]
  security_group_id = aws_security_group.cellenium_study_import_lambda_security_group.id
}

resource "aws_security_group_rule" "cellenium_study_import_security_group_rule_ec2" {
  type                     = "ingress"
  from_port                = var.ec2_security_group_port
  to_port                  = var.ec2_security_group_port
  protocol                 = "TCP"
  source_security_group_id = aws_security_group.cellenium_study_import_lambda_security_group.id
  security_group_id        = var.ec2_security_group_id
  description              = "${data.aws_caller_identity.current.account_id}-${var.stage}-cellenium-lambda-access"
}


resource "aws_lambda_layer_version" "lambda_layer_psycopg2_sqlalchemy" {
  filename   = "${path.module}/lambda_layer.zip"
  layer_name = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.lambda_layer_name}"

  compatible_runtimes      = ["python3.10"]
  compatible_architectures = ["x86_64"]
  depends_on               = [null_resource.lambda_layer_zip]
  source_code_hash         = filebase64sha256("${path.module}/lambda_layer/requirements.txt")
}


resource "aws_lambda_function" "submit_study_import_lambda" {
  filename      = "${path.module}/lambda_submit.zip"
  function_name = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.submit_study_import_lambda_function_name}"
  role          = aws_iam_role.lambda_role.arn
  handler       = "lambda_function.lambda_handler"

  runtime          = "python3.10"
  depends_on       = [data.archive_file.lambda_submit_zip, aws_lambda_layer_version.lambda_layer_psycopg2_sqlalchemy]
  source_code_hash = data.archive_file.lambda_submit_zip.output_base64sha256
  layers           = [aws_lambda_layer_version.lambda_layer_psycopg2_sqlalchemy.arn]

  timeout     = 60
  memory_size = 512

  vpc_config {
    subnet_ids         = var.subnet_ids
    security_group_ids = [aws_security_group.cellenium_study_import_lambda_security_group.id]
  }

  environment {
    variables = {
      AWS_DB_SECRET        = var.db_secret_arn
      BATCH_JOB_DEFINITION = var.batch_job_definition_name
      BATCH_JOB_QUEUE      = var.batch_queue_name
    }
  }
}


resource "aws_lambda_function" "failed_study_import_lambda" {
  filename      = "${path.module}/lambda_failed.zip"
  function_name = "${data.aws_caller_identity.current.account_id}-${var.stage}-${var.failed_study_import_lambda_function_name}"
  role          = aws_iam_role.lambda_role.arn
  handler       = "lambda_function.lambda_handler"

  runtime          = "python3.10"
  depends_on       = [data.archive_file.lambda_failed_zip, aws_lambda_layer_version.lambda_layer_psycopg2_sqlalchemy]
  source_code_hash = data.archive_file.lambda_failed_zip.output_base64sha256
  layers           = [aws_lambda_layer_version.lambda_layer_psycopg2_sqlalchemy.arn]

  timeout     = 60
  memory_size = 512

  vpc_config {
    subnet_ids         = var.subnet_ids
    security_group_ids = [aws_security_group.cellenium_study_import_lambda_security_group.id]
  }

  environment {
    variables = {
      AWS_DB_SECRET = var.db_secret_id
    }
  }
}


resource "aws_lambda_permission" "allow_bucket" {
  statement_id  = "AllowExecutionFromS3Bucket"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.submit_study_import_lambda.arn
  principal     = "s3.amazonaws.com"
  source_arn    = var.s3_bucket_arn
}

resource "aws_s3_bucket_notification" "bucket_notification" {
  bucket = var.s3_bucket_id

  lambda_function {
    lambda_function_arn = aws_lambda_function.submit_study_import_lambda.arn
    events              = ["s3:ObjectCreated:*"]
  }

  depends_on = [aws_lambda_function.submit_study_import_lambda, aws_lambda_permission.allow_bucket]
}


resource "aws_lambda_permission" "allow_cloudwatch_to_call_failed_study_import_lambda" {
  statement_id  = "AllowExecutionFromCloudWatch"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.failed_study_import_lambda.function_name
  principal     = "events.amazonaws.com"
  source_arn    = var.failed_batch_cloudwatch_rule_arn

  depends_on = [aws_lambda_function.failed_study_import_lambda]
}