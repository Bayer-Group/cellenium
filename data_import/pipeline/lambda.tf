resource "null_resource" "lambda_layer_zip" {
  provisioner "local-exec" {
    command = "make build_lambda_layer"
  }
}

resource "null_resource" "lambda_submit_zip" {
  provisioner "local-exec" {
    command = "make build_lambda_submit"
  }
}


resource "null_resource" "lambda_failed_zip" {
  provisioner "local-exec" {
    command = "make build_lambda_failed"
  }
}

resource "aws_security_group" "cellenium_study_import_lambda_security_group" {
  name        = "cellenium_study_import_lambda_security_group"
  description = "Security group for cellenium study import with AWS Batch"
  vpc_id      = data.aws_vpc.vpc.id

  egress {
    from_port       = var.cellenium_ec2_access_security_group_port
    to_port         = var.cellenium_ec2_access_security_group_port
    protocol        = "TCP"
    security_groups = [var.cellenium_ec2_access_security_group_id]
  }
}

resource "aws_security_group_rule" "cellenium_study_import_lambda_security_group_rule" {
  type              = "egress"
  from_port         = 0
  to_port           = 65535
  protocol          = "TCP"
  cidr_blocks       = [data.aws_vpc.vpc.cidr_block]
  security_group_id = aws_security_group.cellenium_study_import_lambda_security_group.id
}

resource "aws_lambda_layer_version" "lambda_layer_psycopg2_sqlalchemy" {
  filename   = "./lambda_layer.zip"
  layer_name = "Python310Psycopg2SQLAlchemy"

  compatible_runtimes      = ["python3.10"]
  compatible_architectures = ["x86_64"]
  depends_on               = [null_resource.lambda_layer_zip]
  source_code_hash         = "$filebase64sha256(./lambda_layer/requirements.zip)"
}

resource "aws_iam_role" "lambda_role" {
  name = "cellenium_study_import_lambda_role"

  assume_role_policy = jsonencode({
    Version   = "2012-10-17"
    Statement = [
      {
        Action    = "sts:AssumeRole"
        Effect    = "Allow"
        Principal = {
          Service = "lambda.amazonaws.com"
        }
      }
    ]
  })

  inline_policy {
    name = "cellenium_study_import_lambda_policy"

    policy = jsonencode({
      Version   = "2012-10-17"
      Statement = [
        {
          Action = [
            "s3:PutObject",
            "s3:GetObject",
            "s3:DeleteObject",
            "s3:ListBucket"
          ]
          Effect   = "Allow"
          Resource = [
            aws_s3_bucket.cellenium_s3_bucket.arn,
            "${aws_s3_bucket.cellenium_s3_bucket.arn}/*"
          ]
        },
        {
          Action = [
            "logs:CreateLogGroup",
            "logs:CreateLogStream",
            "logs:PutLogEvents"
          ]
          Effect   = "Allow"
          Resource = "*"
        },
        {
          Action = [
            "batch:SubmitJob"
          ]
          Effect   = "Allow"
          Resource = "*"
          #          "arn:aws:batch:eu-central-1:850928066808:job-definition/cellenium_data_import:*",
          #                "arn:aws:batch:eu-central-1:850928066808:job-queue/cellenium_data_import_queue"
        },
        {
          Action = [
            "ecr:GetAuthorizationToken",
            "ecr:BatchCheckLayerAvailability",
            "ecr:GetDownloadUrlForLayer",
            "ecr:BatchGetImage"
          ]
          Effect   = "Allow"
          Resource = "*"
        },
        {
          "Effect" : "Allow",
          "Action" : [
            "secretsmanager:GetSecretValue",
            "secretsmanager:DescribeSecret"
          ],
          "Resource" : "arn:aws:secretsmanager:eu-central-1:850928066808:secret:cellenium_database-Rn2qPL"
        }
      ]
    })
  }
}


resource "aws_iam_role_policy_attachment" "lambda_role_policy_attachment_vpc_access" {
  role       = aws_iam_role.lambda_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSLambdaVPCAccessExecutionRole"
}

resource "aws_iam_role_policy_attachment" "lambda_role_policy_attachment_batch_submit" {
  role       = aws_iam_role.lambda_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSBatchServiceEventTargetRole"
}


resource "aws_lambda_function" "submit_study_import_lambda" {
  filename      = "./lambda_submit.zip"
  function_name = "cellenium_submit_study_import"
  role          = aws_iam_role.lambda_role.arn
  handler       = "lambda_function.lambda_handler"

  runtime          = "python3.10"
  depends_on       = [null_resource.lambda_submit_zip, aws_lambda_layer_version.lambda_layer_psycopg2_sqlalchemy]
  source_code_hash = "$filebase64sha256(./lambda_submit/lambda_function.py)"
  layers           = [aws_lambda_layer_version.lambda_layer_psycopg2_sqlalchemy.arn]

  timeout     = 60
  memory_size = 512

    vpc_config {
      # Every subnet should be able to reach an EFS mount target in the same Availability Zone. Cross-AZ mounts are not permitted.
      subnet_ids         = data.aws_subnets.default_vpc_subnets.ids
      security_group_ids = [aws_security_group.cellenium_study_import_lambda_security_group.id]
    }

  environment {
    variables = {
      AWS_DB_SECRET        = var.cellenium_db_secret_name
      BATCH_JOB_DEFINITION = var.cellenium_batch_job_definition_name
      BATCH_JOB_QUEUE      = var.cellenium_batch_job_queue_name
    }
  }
}



resource "aws_lambda_function" "failed_study_import_lambda" {
  filename      = "./lambda_failed.zip"
  function_name = "cellenium_failed_study_import"
  role          = aws_iam_role.lambda_role.arn
  handler       = "lambda_function.lambda_handler"

  runtime          = "python3.10"
  depends_on       = [null_resource.lambda_failed_zip, aws_lambda_layer_version.lambda_layer_psycopg2_sqlalchemy]
  source_code_hash = "$filebase64sha256(./lambda_failed/lambda_function.py)"
  layers           = [aws_lambda_layer_version.lambda_layer_psycopg2_sqlalchemy.arn]

  timeout     = 60
  memory_size = 512

    vpc_config {
      # Every subnet should be able to reach an EFS mount target in the same Availability Zone. Cross-AZ mounts are not permitted.
      subnet_ids         = data.aws_subnets.default_vpc_subnets.ids
      security_group_ids = [aws_security_group.cellenium_study_import_lambda_security_group.id]
    }

  environment {
    variables = {
      AWS_DB_SECRET        = var.cellenium_db_secret_name
    }
  }
}


resource "aws_lambda_permission" "allow_bucket" {
  statement_id  = "AllowExecutionFromS3Bucket"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.submit_study_import_lambda.arn
  principal     = "s3.amazonaws.com"
  source_arn    = aws_s3_bucket.cellenium_s3_bucket.arn
}

resource "aws_s3_bucket_notification" "bucket_notification_h5ad" {
  bucket = aws_s3_bucket.cellenium_s3_bucket.id

  lambda_function {
    lambda_function_arn = aws_lambda_function.submit_study_import_lambda.arn
    events              = ["s3:ObjectCreated:*"]
    filter_suffix       = ".h5ad"
  }

  depends_on = [aws_lambda_function.submit_study_import_lambda, aws_lambda_permission.allow_bucket]
}

resource "aws_s3_bucket_notification" "bucket_notification_h5mu" {
  bucket = aws_s3_bucket.cellenium_s3_bucket.id

  lambda_function {
    lambda_function_arn = aws_lambda_function.submit_study_import_lambda.arn
    events              = ["s3:ObjectCreated:*"]
    filter_suffix       = ".h5mu"
  }

  depends_on = [aws_lambda_function.submit_study_import_lambda, aws_lambda_permission.allow_bucket]
}


resource "aws_lambda_permission" "allow_cloudwatch_to_call_failed_study_import_lambda" {
    statement_id = "AllowExecutionFromCloudWatch"
    action = "lambda:InvokeFunction"
    function_name = aws_lambda_function.failed_study_import_lambda.function_name
    principal = "events.amazonaws.com"
    source_arn = aws_cloudwatch_event_rule.aws_cloudwatch_event_rule.arn

  depends_on = [aws_lambda_function.failed_study_import_lambda]
}